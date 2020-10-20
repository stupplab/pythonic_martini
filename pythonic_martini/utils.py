

import numpy as np
from . import quaternion




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
        
    

    




    