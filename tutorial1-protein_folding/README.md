# Protein Folding using CSP

MELD can fold small proteins using only sequence information (sequence.dat) and secondary structure predictions (ss.dat). 

## Pre-MELD

Before performing the MD simulations, it is recommended to first minimize the system either using Monte Carlo minimizer inside MELD or another MD engine. Here, we used AmberTools to both generate the extended structure from the sequence, and minimize it (tleap.py).

## Set-up

The MELD files will be generated using the setup.py file. To discuss the script, we can divide it into the following sections:

1. MD parameters
```
N_REPLICAS = 30
N_STEPS = 20000
BLOCK_SIZE = 50
```
This block sets the parameters for the REMD simulations. 
