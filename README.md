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



## Development - Adding new code
Any changes/additions you make to the code, do it in a new branch. 
To see which branch you are on, type 
```bash
$ git branch
```
Create a new branch
```bash
$ git branch branchname
```
Go into that branch
```bash
$ git checkout branchname
```
Now this is your branch and you can do any modififications and additions to code. These changes are still in your local workstation. You need to add these to Github repository using `git`. Let's say you added `file1` and modified file `file2.py`. Then do,
```bash
$ git add file1 file2
$ git commit -m "Message telling what these changes are for."
$ git push origin branchname
```
The last command pushes the code to the branch `branchname` in the `origin`, which is already set to be the Github address of the this repository `pythonic_martini`. `branchname` is created if it doesn't already exists.

You can make as make changes as you want using this process but these changes are still in the branch you created and not in the primary branch of this repository called `master`. You can keep using this modified version of the code for your personal use. When and if you think these changes should be incorporated in the master branch, create a **Pull request**. For this, go to the Github page of your branch and click Pull request. Fill up the details asked and submit. On the first time, command `git push origin branchname` also generates a link that you can just copy on your browser to access the pull request page. It may look something like `https://github.com/stupplab/pythonic_martini/pull/new/branchname`.

The maintainer of the repository will now look at the pull request, resolve any conflicts and merge the code as required.

