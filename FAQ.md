# Safe haven for MARTINI/GROMACS questions and tips
Feel free to contribute to this document.

### Tips to reduce trajectory file size
- Dump the trajectory in the form of the compressed xtc file. xtc trajectory does not included velocities and forces, which are usually not needed anyway.
To dump in xtc file, change mdp settings. 
To avoid trr, use 
```
nstxout         = 0
nstvout         = 0
nstfout         = 0
```
To dump xtc,
```
; output frequency
nstxtcout       = 40000 
; xtc precision. 100 implies that coordinates (in nm) are saved upto 2 decimal 
xtc_precision   = 100 
```

- Delete the minimization trajectory file. For instance, min.trr or min.xtc 
are a lot smaller than equilibration trr or xtc but can still occupy resonable space in a high throughout setup consisting hundreds of simulations. min.trr/xtc are not needed generally.


### Where is the equilibrated box dimensions info?
Usually, equilibrated box dimensions will be dumped in the gro file after the equilibration is completed. Many times the simulated termines before completion without dumping the gro file. Since box dimension information in the trajectory file, software `mdtraj` can be used to get the frame wise box dimensions. 
```
import mdtraj
traj = mdtraj.load('PA_water_eq.trr', top='PA_water_min.gro')
Lx, Ly, Lz = traj.unitcell_lengths[-1]
```
There may be other software package that can extract this information