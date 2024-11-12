# Protein Folding using CSP

MELD can fold small proteins using only sequence information (_sequence.dat_) and secondary structure predictions (_ss.dat_). 

## Pre-MELD

Before performing the MD simulations, it is recommended to first minimize the system either using Monte Carlo minimizer inside MELD or another MD engine. Here, we used AmberTools to both generate the extended structure from the sequence, and minimize it (_tleap.py_).

## Set-up

The MELD files will be generated using the _setup.py_ file. To discuss the script, we can divide it into the following sections:

1. MD parameters
```
N_REPLICAS = 30
N_STEPS = 20000
BLOCK_SIZE = 50
```
This block sets the parameters for the REMD simulations. _N_REPLICAS_ is the number of replica ladders, _N_STEPS_ is the total step of simulation, and _BLOCK_SIZE_ dictates how often MELD would save the MD trajectories.
