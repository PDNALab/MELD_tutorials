#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import meld
from meld.remd import ladder, adaptor, leader
import meld.system.montecarlo as mc
from meld import system
from meld.system import patchers
from meld import comm, vault
from meld import parse
from meld import remd
from meld.system import param_sampling
from openmm import unit as u
import mdtraj as md
import glob as glob

N_REPLICAS = 30
N_STEPS = 20000
BLOCK_SIZE = 50

def gen_state(s, index):
    state = s.get_state_template()
    state.alpha = index / (N_REPLICAS - 1.0)
    return state


def get_dist_restraints(filename, s, scaler, ramp, seq):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]
            dist = float(cols[4])
            #
            # note: manually over riding distances
            #
            #dist = 0.45
        
            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=0.0*u.nanometer, r2=0.0*u.nanometer, r3=dist*u.nanometer, r4=(dist+0.2)*u.nanometer, 
                                                 k=350*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:]))
            rest_group.append(rest)
    return dists

def get_dist_restraints_protein(filename, s, scaler, ramp, seq):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]
            dist = float(cols[4])
            #
            # note: manually over riding distances
            #
            #dist = 0.45

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=0.0*u.nanometer, r2=(dist-0.1)*u.nanometer, r3=dist*u.nanometer, r4=(dist+0.1)*u.nanometer,
                                                 k=350*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:]))
            rest_group.append(rest)

    return dists

def setup_system():
    # load the sequence
    sequence = parse.get_sequence_from_AA1(filename='sequence.dat')
    n_res = len(sequence.split())

    # build the system
    p = meld.AmberSubSystemFromPdbFile('protein_min.pdb')
    build_options = meld.AmberOptions(
      forcefield="ff14sbside",
      implicit_solvent_model = 'gbNeck2',
      use_big_timestep = True,
      cutoff = 1.8*u.nanometers,
      remove_com = False,
      #use_amap = False,
      enable_amap = False,
      amap_beta_bias = 1.0,
    )


    builder = meld.AmberSystemBuilder(build_options)
    s = builder.build_system([p]).finalize()
    s.temperature_scaler = system.temperature.GeometricTemperatureScaler(0, 0.3, 300.*u.kelvin, 550.*u.kelvin)

#################

    ramp = s.restraints.create_scaler('nonlinear_ramp', start_time=1, end_time=200,
                                      start_weight=1e-3, end_weight=1, factor=4.0)
    seq = sequence.split()
    for i in range(len(seq)):
        if seq[i][-3:] =='HIE': seq[i]='HIS'
    print(seq)

    #
    # Setup Scaler
    #
    scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    prot_scaler = s.restraints.create_scaler('constant')
    #
    # Torsion Restraints
    #

    #
    # Distance Restraints
    #
    dists = get_dist_restraints('prot_pep_all_fixed.dat', s, scaler, ramp, seq)
    s.restraints.add_selectively_active_collection(dists, int(len(dists)*0.04))
   
    # Keep proteins folded
    prot_rest = get_dist_restraints_protein('ET_contacts.dat',s,prot_scaler,ramp,seq)
    s.restraints.add_selectively_active_collection(prot_rest, int(len(prot_rest)*0.90))

    miniprot_rest = get_dist_restraints_protein('miniprot_contacts.dat',s,prot_scaler,ramp,seq)
    s.restraints.add_selectively_active_collection(miniprot_rest, int(len(prot_rest)*0.90))

    # create the options
    options = meld.RunOptions(
        timesteps = 14286,
        minimize_steps = 20000,
        #min_mc = sched,
        param_mcmc_steps=200
    )

    # create a store
    store = vault.DataStore(s.get_state_template(),N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1, min_acc_prob=0.02)

    remd_runner = remd.leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS,
                                                            ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS, timeout=60000)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms


setup_system()
