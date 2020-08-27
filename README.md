# pythonic_martini

Python scripts for creating and simulating macromolecules using MARTINI / GROMACS, targeted for peptide amphiphiles.

    main.py is the main script.


Original or selectively modified martini files are:\
- martini_v2.0_ions.itp
- martini_v2.2.itp
- martinize.py
- top_all36_prot_forPAs.rtf


Software prerequisites:
- numpy
- vmd 
- GROMACS

Typical simulation steps when using MARTINI to simulation a peptide amphiphile molecule are:
- Make all-atom .pdb file (here using vmd) from a starting pdb file containing generic structure for the molecule
- Prepare coarge-grained files (itp, top, pdb) using martinize.py that uses the all-atom file. 
- Create the simulation box using gromacs filled with coarse-grained molecules
- Solvate the box
- Add ions as needed 
- Run minimization (tiny simulation) to relax the initialization artifacts
- Run equilibration (main simulation)


## Additional details
variable ending with 'filename' are file names with extensions, whereas variables ending with 'name' are without extensions and will be used to read/write from multiple extensions.

