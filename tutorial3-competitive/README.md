# Competitive simulations of two binders

If the binding sites of two binders are similar, we can use MELD to obtain their relative free energies. Since MELD preserves the Boltzmann distribution of structures that are compatible with the data, we can use the relative population of bound binders to calculate ΔΔG.

## Set-up

Tutorial 3 builds on the work done in Tutorial 2. For this stage, we can use the previous python scripts to obtain the following files:

```
1. ET_contacts.dat and miniprot_contacts.dat -- Distances between the CA atoms obtained using pdb_contact.py
2. prot_pep_all_fixed.dat -- N-O distances for ET-miniprotein and ET-TP complexes
3. protein_min.pdb -- Minimized and shifted ET + miniprotein + TP where both the binders are at some distance from ET
```

Similar to Tutorial 2, we apply restraints to keep the two proteins, ET and the miniprotein, folded even at higher replicas. However, unlike in Tutorial 2, the binding modes are already known. Since the binding site is defined by hydrogen bonds between beta strands, we can use N-O distances as the guiding restraints instead of CA,CB-CA,CB pairs. Therefore, for the protein-protein distances, we define three groups:

```
40 N 136 O
40 O 136 N
40 N 88 O
40 O 88 N

42 N 134 O
42 O 134 N
42 N 86 O
42 O 86 N

44 N 132 O
44 O 132 N
44 N 84 O
44 O 84 N
```

In each group, ET (the target) only needs to satisfy one of the restraints essentially allowing it to choose between the two binders. To allow the complex to explore other binding modes, we enforce restraints on just two of these groups.

```python
    dists = get_dist_restraints_strand_pair('prot_pep_all_fixed.dat', s, dist_scaler, ramp, seq)
    s.restraints.add_selectively_active_collection(dists, 2) # Trust two groups
```

## MELD

To prevent co-binding, we introduce three more collections. 

### Bind to only one binder

```python
def exclude_restraint(g,s,i,j,name_i,name_j,seq,scaler=None):
    g.append(s.restraints.create_restraint('distance', scaler,ramp=LinearRamp(0,100,0,1), 
        r1=4.0*u.nanometer, r2=5.0*u.nanometer, r3=7.0*u.nanometer, r4=8.0*u.nanometer, k=250.0*u.kilojoule_per_mole/u.nanometer **2, 
        atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
        atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:])))
```

This collection ensures that only one of the two binders is present in the binding site. Starting from replica 18 and moving down the ladder, we allow ET to bind to one of the binders as dictated by the restraints in the _prot_pep_all_fixed.dat_ file. As we enforce these restraints, the binder that is not chosen should be kept away from the complex. 

We apply two flat-bottom restraints, which are flat between 50 and 70 Å. These restraints are defined between the center of the binding regions of ET and the binder (miniprotein and TP). Since MELD cannot satisfy both the binding and exclusion restraints simultaneously, the unchosen binder must be excluded.

```python
    dist_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.0, alpha_max=0.6, factor=4.0)

    #
    # keep 1 of TP or miniprotein far from ET
    #
    exclude_all = []

    # TP
    exclude_rests = []
    exclude_restraint(exclude_rests,s,ET_center[0],binding_centers[0],'CA','CA',seq,scaler=dist_scaler)
    exclude_all.append(s.restraints.create_restraint_group(exclude_rests,1))

    # miniprotein
    exclude_rests = []
    exclude_restraint(exclude_rests,s,ET_center[0],binding_centers[1],'CA','CA',seq,scaler=dist_scaler)
    exclude_all.append(s.restraints.create_restraint_group(exclude_rests,1))

    s.restraints.add_selectively_active_collection(exclude_all, 1)
```

### Prevent co-binding

```python
def non_interacting(g,s,i,j,name_i,name_j,seq,scaler=None):
    g.append(s.restraints.create_restraint('distance', scaler,ramp=LinearRamp(0,100,0,1), 
        r1=0.0*u.nanometer, r2=4.5*u.nanometer, r3=100.0*u.nanometer, r4=101.0*u.nanometer, k=250.0*u.kilojoule_per_mole/u.nanometer **2,
        atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
        atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:])))
```

One way to satisfy both the binding and exclusion restraints is by forming a megaplex where all three components are bound to each other. To prevent this, we can add another restraint that is flat only when the two binders are within a distance of 45 to 100 Å from each other.

```python
    exclude_rests = []
    non_interacting(exclude_rests,s,pep_centers[0],pep_centers[1],'CA','CA',seq,scaler=const_scaler)
    s.restraints.add_as_always_active_list(exclude_rests)
```

### Keep excluded binder within some distance from ET

Finally, we add restraints to ensure that both binders remain within a certain distance from ET and prevent it from drifting too far away.

```python
    scaler3 = s.restraints.create_scaler('constant')
    conf_rest = []
    for i in binding_centers:
        conf_rest.append(s.restraints.create_restraint('distance', scaler3,ramp=LinearRamp(0,100,0,1),
            r1=0.0*u.nanometer, r2=0.0*u.nanometer, r3=7.0*u.nanometer, r4=8.0*u.nanometer, k=250.0*u.kilojoule_per_mole/u.nanometer **2,
            atom1=s.index.atom(i,'CA', expected_resname=seq[i][-3:]),
            atom2=s.index.atom(ET_center[0],'CA', expected_resname=seq[ET_center[0]][-3:])))

    s.restraints.add_as_always_active_list(conf_rest)
```

## Analysis

Similar to Tutorials 1 and 2, we generate the MELD files by executing _setup.py_. A simple way to analyze the trajectories is to calculate N-O distances between ET and the two binders at each replica.

```
parm protein.prmtop
reference protein_min.pdb

trajin ../trajectories/trajectory.00.dcd 9999 19998

rms reference out ET_rms.dat :1-68@N,O

# miniprot
distance mini1 :145@N,O :44@N,O out mini_ET_1.00.dat
distance mini2 :147@N,O :42@N,O out mini_ET_2.00.dat
distance mini3 :149@N,O :40@N,O out mini_ET_3.00.dat

# tp
distance tp1 :84@N,O :44@N,O out TP_ET_1.00.dat
distance tp2 :86@N,O :42@N,O out TP_ET_2.00.dat
distance tp3 :88@N,O :40@N,O out TP_ET_3.00.dat

go
```

We can then calculate the average of these distances at each timestep and set a threshold (in this case, 4 Å) to determine if the binder is bound. By counting how often this threshold is met for both binders, we can calculate their relative populations in the bound state.

```
for((i=0;i<=9;i++)); do
    #echo cpptraj -i distance.in
    sed "s/00/0$i/g" distance.in > distance_$i.in
    cpptraj -i distance_$i.in
    paste -d '=' mini_ET_{1..3}.0$((i)).dat | sed 's/=\S*//g' | awk 'NR>1 {print ($2 + $4 +$6)/3}' > mini_dist_ave.0$((i)).dat
    paste -d '=' TP_ET_{1..3}.0$((i)).dat | sed 's/=\S*//g' | awk 'NR>1 {print ($2 + $4 +$6)/3}' > TP_dist_ave.0$((i)).dat
    awk '$0 < 4' mini_dist_ave.0$((i)).dat | wc -l | awk '{print $0/10000}' >> mini_dist_percent.4.dat
    awk '$0 < 4' TP_dist_ave.0$((i)).dat | wc -l | awk '{print $0/10000}' >> TP_dist_percent.4.dat
done

for((i=10;i<=29;i++)); do
    #echo cpptraj -i distance.in
    sed "s/00/$i/g" distance.in > distance_$i.in
    cpptraj -i distance_$i.in
    paste -d '=' mini_ET_{1..3}.$((i)).dat | sed 's/=\S*//g' | awk 'NR>1 {print ($2 + $4 +$6)/3}' > mini_dist_ave.$((i)).dat
    paste -d '=' TP_ET_{1..3}.$((i)).dat | sed 's/=\S*//g' | awk 'NR>1 {print ($2 + $4 +$6)/3}' > TP_dist_ave.$((i)).dat
    awk '$0 < 4' mini_dist_ave.$((i)).dat | wc -l | awk '{print $0/10000}' >> mini_dist_percent.4.dat
    awk '$0 < 4' TP_dist_ave.$((i)).dat | wc -l | awk '{print $0/10000}' >> TP_dist_percent.4.dat
done

rm distance_*.in
```
