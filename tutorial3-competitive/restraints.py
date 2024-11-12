import numpy as np
import meld
from meld.remd import ladder, adaptor, leader
import meld.system.montecarlo as mc
from meld import system
from meld.system import patchers
from meld import comm, vault
from meld import parse, util
from meld import remd
from meld.system import param_sampling
from openmm import unit as u
import mdtraj as md
import glob as glob
from meld.system.scalers import LinearRamp,ConstantRamp

def non_interacting(g,s,i,j,name_i,name_j,seq,scaler=None):
    g.append(s.restraints.create_restraint('distance', scaler,ramp=LinearRamp(0,100,0,1), 
        r1=0.0*u.nanometer, r2=4.5*u.nanometer, r3=100.0*u.nanometer, r4=101.0*u.nanometer, k=250.0*u.kilojoule_per_mole/u.nanometer **2,
        atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
        atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:])))

def exclude_restraint(g,s,i,j,name_i,name_j,seq,scaler=None):
    g.append(s.restraints.create_restraint('distance', scaler,ramp=LinearRamp(0,100,0,1), 
        r1=4.0*u.nanometer, r2=5.0*u.nanometer, r3=7.0*u.nanometer, r4=8.0*u.nanometer, k=250.0*u.kilojoule_per_mole/u.nanometer **2, 
        atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
        atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:])))

# Restraints on Protein CA
def make_cartesian_collections(s, scaler, residues, ramp, delta=0.35*u.nanometer, k=250.*u.kilojoule_per_mole / u.nanometer**2):
    cart = []
    backbone = ['CA']
    #Residues are 1 based
    #index of atoms are 1 base
    for i in residues:
        #print(i)
        for b in backbone:
            #print(b)
            atom_index = s.index.atom(i,b) 
            x,y,z = s.template_coordinates[atom_index]/10.
            x=x*u.nanometer; y=y*u.nanometer; z=z*u.nanometer
            rest = s.restraints.create_restraint('cartesian',scaler, LinearRamp(0,15,0,1), atom_index=atom_index,
                x=x, y=y, z=z, delta=delta,force_const=k)
            cart.append(rest)
    return cart

