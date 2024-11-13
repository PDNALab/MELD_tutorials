# Modelling Protein-Protein complexes using NMR CSP data

MELD can verify the binding modes of protein-protein complexes. Here, we will use NMR CSP data to extract the possible binding residues of the target protein. Residues with the largest CSP values, specifically those with ΔδN,H greater than 0.25 ppm, were considered _active residues_ and are presumed to be near the actual binding site. These values were calculated for ET-peptide complexes but are considered transferable to other ET-protein complexes. More information about these values can be found in:

[Mondal, A.; Swapna, G.; Lopez, M. M.; Klang, L.; Hao, J.; Ma, L.; Roth, M. J.; Montelione, G. T.; Perez, A. Structure Determination of Challenging Protein–Peptide Complexes Combining NMR Chemical Shift Data and Molecular Dynamics Simulations. J. Chem. Inf. Model. 2023, 63, 2058–2072](https://pubs.acs.org/doi/10.1021/acs.jcim.2c01595)

## Pre-MELD

In addition to minimizing the system, we have three additional steps before executing the _setup.py_ file.

### Generate protein-protein pairs

We will use the _process_CSP.py_ script to generate the restraint file for the protein-protein pairs. The _active residues_ are defined in the list named _final_, and the sequences of the target and the binder are defined by the strings _ET_ and _minprot_, respectively.

```python
final=[3,6,9,20,23,24,25,33,35,36,37,38,41,42,43,44,45,46,47]

ET=str('SYDEKRQLSLDINRLPGEKLGRVVHIIQSREPSLRDSNPDEIEIDFETLKPTTLRELERYVKSCLQKK') #TP
minprot=str('MKEEEIEEMIEKAKEELRKRYPEAKEVFLSFTYEVNGKLKIKLTRFDPSMSLEEVEERIEEEVKRLLKEADSIEIRVHTTV')
seq=ET+minprot
```

As with the target protein, we also need to define the potential binding residues of the binder. In this example, we assume the binding mode involves the beta-strand residues, specifically residues 25 to 47.

```python
hairpin1=25
hairpin2=47
```

We can then write our restraint file which is basically all the possible CA,CB-CA,CB combinations between the _active residues_ in the target and in the binder.

```python
file=open('prot_pep_all_fixed.dat','w')
for i in final:
    for j in range(69,150): # add 1
        if j < hairpin1+69 or j > hairpin2+69:
            continue
        if seq[i-1] =='G' and seq[j-1] != 'G':
            file.write( "{} CA {} CB {}\n".format(i, j, 0.8 ))
            file.write("\n")
        if seq[i-1] !='G' and seq[j-1] =='G':
            file.write( "{} CB {} CA {}\n".format(i, j, 0.8 ))
            file.write("\n")
        if seq[i-1] =='G' and seq[j-1] =='G':
            file.write( "{} CA {} CA {}\n".format(i, j, 0.8 ))
            file.write("\n")
        if seq[i-1] !='G' and seq[j-1] !='G':
            file.write( "{} CB {} CB {}\n".format(i, j, 0.8 ))
            file.write("\n")
    
file.close()
```

### Keeping the proteins folded

Since we are only testing the binding modes of the protein-protein complex, it is best to keep the proteins folded even at high replica indices. To achieve this, we can generate distance restraints that maintain the conformations of both the target and the binder. This is done using pdb_contact.py, which generates a file containing all CA-CA pairs within a PDB file along with their respective distances.

### Starting conformation

Our initial conformation should not be biased toward any specific binding mode. To ensure this, the peptide must be at least 30 A away from the receptor. This can be achieved using the change_coor.py script, which shifts the coordinates of all atoms to be 30 A from their positions in the original PDB file.

## Set-up

Discussions of the MELD options can be found in Tutorial 1. Here, we will focus on how the restraints were defined. We have two collections for this simulation -- (1) collections to keep the proteins folded, and (2) collection to drive protein-protein complex formation.

### Keeping the proteins folded

This function reads the _ET_contacs.dat_ and _miniprot_contacts.dat_ files to generate the distance restraints for the proteins. In addition to the atom and residue information, these files also have the distance values as the fifth column. Our restraint is flat between (dist-1) A and dist A, ensuring that the proteins only adopt one conformation.

```python
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
```

We then define two collections, one for each protein, where we satisy 90% of the restraints. These restraints have a constant force constant in all the replicas.

```python
    prot_scaler = s.restraints.create_scaler('constant')

    # Keep proteins folded
    prot_rest = get_dist_restraints_protein('ET_contacts.dat',s,prot_scaler,ramp,seq)
    s.restraints.add_selectively_active_collection(prot_rest, int(len(prot_rest)*0.90))

    miniprot_rest = get_dist_restraints_protein('miniprot_contacts.dat',s,prot_scaler,ramp,seq)
    s.restraints.add_selectively_active_collection(miniprot_rest, int(len(prot_rest)*0.90))
```

### Protein-protein restraints

For the protein-protein restraints, we apply weaker contacts, with restraints flat between 0 and 8 A to allow the proteins to thoroughly sample different binding modes.

```python
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
```

Since we include all possible CA-CB pairs between the proteins, only a fraction of these restraints need to be satisfied to form the protein-protein complex. In this case, only 4% of the restraints in the _prot_pep_all_fixed.dat_ file are required to be satisfied.

```python
    scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)

    dists = get_dist_restraints('prot_pep_all_fixed.dat', s, scaler, ramp, seq)
    s.restraints.add_selectively_active_collection(dists, int(len(dists)*0.04))
```

## Analysis

Similar to Tutorial 1, we will generate the MELD files by executing _setup.py_. After the simulations, we can extract the trajectories using the _extract_trajs.slurm_ file.
