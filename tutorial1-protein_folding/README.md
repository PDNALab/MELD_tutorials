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

2. CPI Restraints

For this example, the folding is guided by two insights: (1) proteins form a hydrophobic core, and (2) beta-strands can form hydrogen bonds with other strands. We also have a third restraint based on the insight that proteins have to be compact to fold. For insight 1, we'll have to first define the hydrophobic residues:
```
hydrophobes = 'AILMFPWV'
hydrophobes_res = ['ALA','ILE','LEU','MET','PHE','PRO','TRP','VAL']
```
Insight 2, on the other hand, depends on the provided secondary structure information. Hydrogen bonds, in this case, is only defined for residues that do not belong in the same strand. From these two information, we can then generate the _hydrophobe.dat_ and _strand_pair.dat_ files which contains the atoms that will form hydrophobic and strand pairings.
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
