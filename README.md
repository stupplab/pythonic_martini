# pythonic_martini

Python scripts for creating and simulating macromolecules using MARTINI / GROMACS, targeted for peptide amphiphiles.

    Built using python3.8


## Software prerequisites:
- Python3 
- Python2 (for martinize.py)
- numpy
- vmd     
- GROMACS

(vmd and gromacs should be callable at the command line)

You can download the package and run the example script using the following commands. Running example_script will create the directory 'example', then prepare and run the PA-peptide coassembly simulation inside it.
    
    git clone https://github.com/stupplab/pythonic_martini.git
    cd pythonic_martini
    python example_script.py


If you want to run the example script from anywhere, i.e., import pythonic_martini from anywhere, install the package in you local python package directory by running

    python setup.py install --user


Typical simulation steps when using MARTINI to simulation a peptide amphiphile molecule are:
- Make all-atom .pdb file (here using vmd) from a starting pdb file containing generic structure for the molecule
- Prepare coarge-grained files (itp, top, pdb) using martinize.py that uses the all-atom file. 
- Create the simulation box using gromacs filled with coarse-grained molecules
- Solvate the box
- Add ions as needed 
- Run minimization (tiny simulation) to relax the initialization artifacts
- Run equilibration (main simulation)

main.py is the main script that lets you do all the above as shown by the example



## Additional details
variable ending with 'filename' are file names with extensions, whereas variables ending with 'name' are without extensions and will be used to read/write from multiple extensions.

Original or selectively modified martini files are:\
- martini_v2.0_ions.itp
- martini_v2.2.itp
- martinize.py
- top_all36_prot_forPAs.rtf


