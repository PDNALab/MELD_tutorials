# Protein Folding using CSP

MELD can fold small proteins using only sequence information (_sequence.dat_) and secondary structure predictions (_ss.dat_). 

## Pre-MELD

Before performing the MD simulations, it is recommended to first minimize the system either using Monte Carlo minimizer inside MELD or another MD engine. Here, we used AmberTools to both generate the extended structure from the sequence, and minimize it (_tleap.py_).

## Set-up

The MELD files will be generated using the _setup.py_ file. To discuss the script, we can divide it into the following sections:

### MD Parameters
   
```
N_REPLICAS = 30
N_STEPS = 20000
BLOCK_SIZE = 50
```

This block sets the parameters for the REMD simulations. _N_REPLICAS_ is the number of replica ladders, _N_STEPS_ is the total step of simulation, and _BLOCK_SIZE_ dictates how often MELD would save the MD trajectories.

### CPI Restraints

For this example, the folding is guided by the following insights: proteins form a hydrophobic core, and beta-strands can form hydrogen bonds with other strands. We also have a third restraint based on the insight that proteins have to be compact to fold. For insight 1, we'll have to first define the hydrophobic residues:

```
hydrophobes = 'AILMFPWV'
hydrophobes_res = ['ALA','ILE','LEU','MET','PHE','PRO','TRP','VAL']
```

Insight 2, on the other hand, depends on the provided secondary structure information. Hydrogen bonds, in this case, is only defined for residues that do not belong in the same strand. From these two information, we can then generate the _hydrophobe.dat_ and _strand_pair.dat_ files which contain the atoms that will form hydrophobic and strand pairings, respectively.

```python
def create_hydrophobes(s,group_1=np.array([]),group_2=np.array([]),CO=True):
    hy_rest=open('hydrophobe.dat','w')
    atoms = {"ALA":['CA','CB'],
             "VAL":['CA','CB','CG1','CG2'],
             "LEU":['CA','CB','CG','CD1','CD2'],
             "ILE":['CA','CB','CG1','CG2','CD1'],
             "PHE":['CA','CB','CG','CD1','CE1','CZ','CE2','CD2'],
             "TRP":['CA','CB','CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
             "MET":['CA','CB','CG','SD','CE'],
             "PRO":['CD','CG','CB','CA']}
    #Groups should be 1 centered
    n_res = s.residue_numbers[-1]
    print(n_res)
    group_1 = group_1 if group_1.size else np.array(list(range(n_res)))+1
    group_2 = group_2 if group_2.size else np.array(list(range(n_res)))+1

    #Get a list of names and residue numbers, if just use names might skip some residues that are two
 #times in a row
    #make list 1 centered
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    print(sequence)
    sequence = dict(sequence)

    print(group_1)
    print(group_2)
    group_1 = [ res for res in group_1 if (sequence[res-1] in hydrophobes_res) ]
    group_2 = [ res for res in group_2 if (sequence[res-1] in hydrophobes_res) ]

    print(group_1)
    print(group_2)
    pairs = []
    hydroph_restraints = []
    for i in group_1:
        for j in group_2:

            # don't put the same pair in more than once
            if ( (i,j) in pairs ) or ( (j,i) in pairs ):
                continue

            if ( i ==j ):
                continue

            if (abs(i-j)< 7):
                continue
            pairs.append( (i,j) )

            atoms_i = atoms[sequence[i-1]]  #atoms_i = atoms[sequence[i]]
            atoms_j = atoms[sequence[j-1]]  #atoms_j = atoms[sequence[j]]

            local_contact = []
            for a_i in atoms_i:
                for a_j in atoms_j:
                    hy_rest.write('{} {} {} {}\n'.format(i,a_i, j, a_j))
            hy_rest.write('\n')

def generate_strand_pairs(s,sse,subset=np.array([]),CO=True):
    f=open('strand_pair.dat','w')
    n_res = s.residue_numbers[-1]
    subset = subset if subset.size else np.array(list(range(n_res)))+1
    strand_pair = []
    for i in range(len(sse)):
        start_i,end_i = sse[i]
        for j in range(i+1,len(sse)):
            start_j,end_j = sse[j]

            for res_i in range(start_i,end_i+1):
                for res_j in range(start_j,end_j+1):
                    if res_i in subset or res_j in subset:
                        f.write('{} {} {} {}\n'.format(res_i, 'N', res_j, 'O'))
                        f.write('{} {} {} {}\n'.format(res_i, 'O', res_j, 'N'))
                        f.write('\n')
```

These files will only generate text files for debugging purposes. To convert these to distance restraints and python blocks, we need the following functions. These define flat-bottom restraints where no energy penalty is added between r2 and r3. Since strand pairing is guided by hydrogen bonds, we enfore stronger contacts for these restraints (i.e., flat until 3.5 A only compared to 5 A for hydrophobes).

```python
def get_dist_restraints_hydrophobe(filename, s, scaler, ramp, seq):
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

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=0.0*u.nanometer, r2=0.0*u.nanometer, r3=0.5*u.nanometer, r4=0.7*u.nanometer,
                                                 k=250*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:]))
            rest_group.append(rest)
    return dists

def get_dist_restraints_strand_pair(filename, s, scaler, ramp, seq):
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
            #

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=0.0*u.nanometer, r2=0.0*u.nanometer, r3=0.35*u.nanometer, r4=0.55*u.nanometer,
                                                 k=250*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:]))
            rest_group.append(rest)
    return dists
```

### Set-up

Once we have set up these functions, we can them call them to generate the MELD files.
```python
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
```
This block generates the AMBER topology for both the protein and the solvent, and defines the timestep (_use_big_timestep_ being enabled as True implies we are using a timestep of 3.5 fs). This also defines the temperature along the replica ladder. We use a geometric scaler for the temperature between replicas 1 to 9 with the range of 300 K to 550 K. For replicas 10 and beyond, the temperature is constant at 550 K.

In contrast to temperature, we only scale the restraints from replicas 12 to 30. We use a nonlinear scaler for the distance which enforces maximum restraints at replica 12 and zero restraints at replica 30. We also define here the accuracy of the collections. For hydrophobic pairing, we only enforce 1.2 * N<sub>h</sub> contacts where N<sub>h</sub> is the number of hydrophobic residues. Strand pairing, on the other hand, only enforces 0.45 * N<sub>e</sub> restraints where N<sub>e</sub> is the number of residues predicted to be a beta-strand.

```python
    #
    # Setup Scaler
    #
    scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    subset1= np.array(list(range(n_res))) + 1
    #


    create_hydrophobes(s,group_1=subset1,group_2=subset1,CO=False)
    
    
    ##hydrophobic contacts
    dists = get_dist_restraints_hydrophobe('hydrophobe.dat', s, scaler, ramp, seq)
    s.restraints.add_selectively_active_collection(dists, int(1.2 * no_hy_res))   

    ##strand pairing
    sse,active = make_ss_groups(subset=subset1)
    generate_strand_pairs(s,sse,subset=subset1,CO=False)
    #
    dists = get_dist_restraints_strand_pair('strand_pair.dat', s, scaler, ramp, seq)
    s.restraints.add_selectively_active_collection(dists, int(0.45*active))
```

### Replica Exchange

The blocks specify some options for REMD. Since we used _big_timestep_, the 14286 value for the timesteps indicates that we are attempting exchanges between the replicas every 50 ps. One exchange attempt corresponds to one REMD step. For this example, we are doing 20000 steps which is equivalent to 1000 ns of MD simulation.

```python
    # create the options
    options = meld.RunOptions(
        timesteps = 14286,
        minimize_steps = 20000,
        min_mc = sched,
        param_mcmc_steps=200
    )
```

Other options for REMD are found in this block. It is in our best interest to allow more exchange attempts during the simulation. The _n_trials_ parameter controls this; here, we perform 48×48 attempts to exchange adjacent replicas.

```python
    # create a store
    store = vault.DataStore(gen_state(s,0), N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)  # why i need gen_state(s,0)? doubtful
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1, min_acc_prob=0.02)

    remd_runner = remd.leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS, timeout=60000)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()
```

## Analysis

The _setup.py_ can be executed to the Data folder which contains the MELD files. To run the simulation, we can submit a SLURM job using the _run.slurm_ file. After the simulation is completed, the trajectories can be generated using the _extract_trajs.slurm_ file.
