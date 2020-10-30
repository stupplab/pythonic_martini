

import numpy as np
from . import quaternion




def actual_atoms_added(filename, keyword):
    """Count the number of lines with <keyword> in the <filename> - usually gro or pdb file
    """
    num_keyword = 0
    with open(filename, 'r') as f:
        for line in f:
            if keyword in line:
                num_keyword += 1

    return num_keyword




def actual_molecules_added(filename, itpfilename, start_linenumber=0):
    """Calculate the number of molecules, defined in <itpfilename>, 
    in <filename> - usually gro or pdb file
    """
    
    # get atom names from the <molname>.itp
    atom_names=[]
    start=False
    with open(itpfilename, 'r') as f:
        for line in f:
            if not start:
                if '[ atoms ]' in line:
                    start = True
                    continue
            if start:
                if line.split()==[]:
                    start = False
                    break
                atom_names += [line.split()[3]]
      
    # get all the atom names from the <outfilename>
    with open(filename, 'r') as f:
        lines = f.readlines()
    all_names = [] # e.g. PAM, ALA, GLU
    for line in lines[start_linenumber:]:
        if line.split()==[]:
            continue
        all_names += [line.split()[0].lstrip('0123456789')]

    # now count number of atom_names sequence in all_names
    actual_nmol = 0
    natoms = len(atom_names)
    i = 0
    while i<len(all_names):
        if all_names[i:i+natoms] == atom_names:
            i += natoms
            actual_nmol += 1
        else:
            i += 1


    return actual_nmol




def make_molecule_template(mol_names, mol_numbers, atom_names, atom_positions, filename):
    '''Create the molecule scaffold structure to be used by VMD to make the all-atom structure
    molecules: list of [molecule name, position list]
    in order to how it should appear in the structure file. 
    '''


    lines = ['COMPND    UNNAMED']
    for i,(mol_name,mol_number,atom_name,atom_position) in enumerate(zip(mol_names,mol_numbers,atom_names,atom_positions)):

        atom_number = i+1
        lines += ['ATOM'+ '{:>7}  '.format(atom_number)+   '{:<4}'.format(atom_name)+ '{}'.format(mol_name)+ '{:>6}  '.format(mol_number)+ '{:>8.2f}'.format(atom_position[0])+ '{:>8.2f}'.format(atom_position[1])+ '{:>8.2f}'.format(atom_position[2])+ '  1.00  0.00'+  '{:>12}'.format(atom_name[0])]

    lines += ['END']
    
    
    with open(filename, 'w') as f:
        f.write('\n'.join(lines))
        
    

    




    