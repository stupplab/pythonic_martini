"""
Script that runs a co-assembly simulation of a peptide amphiphile 
and a peptide using MARTINI-GROMACS
"""
#########################################


import numpy as np
import os
from pythonic_martini import main




## Simulation parameters
PA_seq             = 'C16VVAAEE'
num_PA             = 4
pep_seq            = 'DFFPV'
num_pep            = 8
Ctermini_type      = 'OH'           # currently accepts OH or NH2 
PA_charge          = -1             # Net charge on the PA molecule - places the charge at the outermost possible residue
pep_charge         = -1
Lx                 = 4
Ly                 = 4
Lz                 = 4



################## Now the script that prepares and runs the simulation ##################
## Create a directory for this example and go inside it
if not os.path.exists('example'):
    os.mkdir('example')
os.chdir('example')


## Make the all-atom <name>_aa.pdb of PA and pep using vmd, the generic PA backbone structure used is 
## already created. See PA_generic.pdb
main.make_aa_pdb(PA_seq, name='PA')
main.make_aa_pdb(pep_seq, name='pep')


## Create coarse-grained version of PA and pep molecules using the <name>_aa.pdb
## This generate <name>.itp, <name>.top and <name>.pdb files
main.create_CGfiles_using_martinizepy(Ctermini_type, PA_charge, name='PA')
main.create_CGfiles_using_martinizepy(Ctermini_type, pep_charge, name='pep')


## Create a simulation box with desired number of PA molecules
## Put the molecule info in a new .top file
main.create_simulation_box(Lx,Ly,Lz, 
    nmol        = num_PA,
    molname     = 'PA',
    vdwradius   = 0.4,
    boxfilename = 'coPA_box.gro',     
    topfilename = 'coPA.top')


## Insert pep molecules in the box and update the .top file
main.insert_molecules('coPA_box.gro', 'pep', 
    nmol        = num_pep, 
    vdwradius   = 0.4, 
    outfilename = 'coPA_box.gro', 
    topfilename = 'coPA.top')


## Solvate the simulation box
main.insert_water('coPA_box.gro', 
    volume      = Lx*Ly*Lz,
    vdwradius   = 0.15,
    outfilename = 'coPA_water.gro',
    topfilename = 'coPA.top')



"""# another function to solvate the simulation box
main.solvate('coPA_box.gro', 
    vdwradius   = 0.21,
    topfilename = 'coPA.top',
    outfilename = 'coPA_water.gro',
    )
"""


## Now neutralize the charge in the system due to PA and pep
## Add NA+ or CL- ions depending on the charge

main.neutralize_system(
    inwhichfilename = 'coPA_water.gro',
    molname         = 'PA',
    vdwradius       = 0.21,
    outfilename     = 'coPA_water.gro',
    topfilename     = 'coPA.top')

main.neutralize_system(
    inwhichfilename = 'coPA_water.gro',
    molname         = 'pep',
    vdwradius       = 0.21,
    outfilename     = 'coPA_water.gro',
    topfilename     = 'coPA.top')



## Increation the ionic strength (if desired) by adding <nions> amount of NA and CL atoms.

main.add_ions(
    inwhichfilename = 'coPA_water.gro',
    nions           = 10, 
    vdwradius       = 0.21, 
    outfilename     = 'coPA_water.gro',
    topfilename     = 'coPA.top')


## copy the minimize.mdp (can be custom made) file into the example directory
os.system('cp %s/minimize.mdp ./'%os.path.dirname(main.__file__))

## Makes the tpr file and runs the minimization of the simulation
main.minimization(
    mdpfilename    = 'minimize.mdp',
    topfilename    = 'coPA.top',
    grofilename    = 'coPA_water.gro',
    tprname        = 'coPA_water_min')


## Copy the equilibrate.mdp (can be custom made) file into the example directory
os.system('cp %s/equilibrate.mdp ./'%os.path.dirname(main.__file__))

## Prepare the .tpr file to be used by the gromacs mdrun
main.prepare_eq_tpr(
    mdpfilename    = 'equilibrate.mdp',
    topfilename    = 'coPA.top',
    grofilename    = 'coPA_water_min.gro',
    tprname        = 'coPA_water_eq')

## Run the simulation
main.equilibration('coPA_water_eq')


## Remove any backup files created by the gromacs in this whole process
os.system('rm \#*')

