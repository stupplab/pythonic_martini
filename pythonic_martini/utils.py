

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
        
    

    

def init_fiber_config_co(grofile, num_atomss, num_moleculess, Lx, Ly, Lz, start_from_nth_atoms, inverts, BB_indicess,
    delta = 0.5, radial_offset = 0.5):
    """ Change the random positioning in grofile with fiber like initial config
    assumes all molecules are in contiguous in the grofile
    Note: start_from_nth_atom starts from index 0, whereas in grofile atom index starts from 1.
    BB is backbone.
    Assumes first residue is BB
    BB_indices: Carbon indices to use when aligning. Indices are after invertion if invert=True
    Vector BB_indices[1]-BB_indices[0] is used.
    If invert is true, the other end is pointed towards fiber inside

    num_atoms, mol_start_atom_indices, invert, C_indices are lists/arrays corresponding to grofiles' molecules
    which_mol are indices of molecules in grofiles corresponding to mol_start_atom_indices
    mol_indices are molecule number in the grofile

    delta = 0.5 # distance between layers in nm
    radial_offset = 0.5 # in nm

    NOTE: VERY ARDOUS CODE. NOT HAPPY
    """


    # Parameters for the initial configuration as a fiber 
    num_PA_layer = 9
    theta = np.pi*40/180 # angle between two PA in a layer
    theta_offset = np.pi*20/180
    

    atoms_positionss = []
    BB_idss = []  # C_ids contain atom number of C in atom positions
    for k,num_atoms in enumerate(num_atomss):
        # Molecule's atoms' positions. Read from grofile
        atoms_positions = []
        invert = inverts[k]
        BB_indices = BB_indicess[k]
        start_from_nth_atom = start_from_nth_atoms[k]
        BB_ids = []
        with open(grofile, 'r') as f:
            lines = f.readlines()
            start = False
            BB_pos = []
            for line in lines[2:]:
                words = line.split()
                try:
                    float(words[-1])
                except:
                    continue
                if len(words)<4:
                    continue

                index = int(line[-30:-25])
                
                if index == start_from_nth_atom+1:
                    start = True
                if start:
                    if len(atoms_positions) == num_atoms:
                        break
                    
                    atoms_positions += [ [float(words[-3]), float(words[-2]), float(words[-1])] ]
                    if 'BB' == words[1][0:2]:
                        BB_pos += [ atoms_positions[-1] ]
                        BB_ids += [len(atoms_positions)-1]
        if invert:
            atoms_positions = np.array(atoms_positions)#[::-1]
            BB_ids = BB_ids[::-1]
        else:
            atoms_positions = np.array(atoms_positions)
        

        # Translate with origin at outermost (C16) Carbon of alkyl tail
        atoms_positions -= atoms_positions[BB_ids[BB_indices[0]]]

        
        # Align atoms_positions in the x-y plane by aligning the vector
        # between first two consecutive C atoms of C16 towards x-axis
        # NOTE: Assumes the first atoms is C16
        v = np.array(atoms_positions[BB_ids[BB_indices[1]]]) - np.array(atoms_positions[BB_ids[BB_indices[0]]])
        q = quaternion.q_between_vectors(v, [1,0,0])
        for i,pos in enumerate(atoms_positions):
            atoms_positions[i] = quaternion.qv_mult(q, pos)
            
        atoms_positionss += [atoms_positions]
        BB_idss += [BB_ids] 
    
        
    BB_positions = []  # of the outermost C  or innermost (of the core) C
    quaternions = [] # wrt when PAM chain is oriented in x-axis
    for i in range(sum(num_moleculess)):
        j = i  % num_PA_layer
        k = i // num_PA_layer
        
        q = quaternion.axisangle_to_q([0,0,1], theta * j + theta_offset * (1+(-1)**k)/2 )
        quaternions += [q]

        pos = np.array([0,0,k*delta]) + quaternion.qv_mult(q, [radial_offset,0,0])
        BB_positions += [ pos ]
    

    # randomly arrange molecules
    args = np.arange(sum(num_moleculess))
    np.random.shuffle(args)
    bins = []
    for arg in args:
        s = 0
        for i,num_atoms in enumerate(num_atomss):
            if arg < s + num_moleculess[i]:
                bins+=[i]
                break
            else:
                s += num_moleculess[i]
    

    # Calculate positions of all atoms
    t = np.array(num_atoms)*np.array(num_moleculess)
    positions = []
    for i,arg in enumerate(args):
        atoms_positions = atoms_positionss[bins[i]]
        invert = inverts[bins[i]]
        positions_ = []
        for j,pos in enumerate(atoms_positions):
            p = quaternion.qv_mult(quaternions[i], pos)
            p += BB_positions[i]
            positions_ += [ p ]
        # if invert:
        #     positions += positions_#[::-1]
        # else:
        #     positions += positions_
        positions += [ positions_ ]


    positions_ = [[]]*len(positions)
    for i,pos in enumerate(positions):
        positions_[args[i]] = pos
    positions = positions_
    
    positions_=[]
    for pos in positions:
        positions_ += pos
    positions = np.array(positions_)
    

    positions = np.array(positions) - np.mean(positions, axis=0) +[Lx/2,Ly/2,Lz/2]
    

    # Write new atom positions in the grofile
    new_lines = []
    n=0
    start_from_nth_atom = 0 
    with open(grofile, 'r') as f:
        lines = f.readlines()

        for line in lines:
            words = line.split()
            try:
                float(words[-1])
            except:
                new_lines += [line]
                continue
            if len(words)<4:
                new_lines += [line]
                continue
            
            index = int(line[-30:-25])
            if (index < start_from_nth_atom+1) or (n>=len(positions)):
                new_lines += [line]
            else:
                line_ = line[:-25]+'%8.3f%8.3f%8.3f\n'%tuple(positions[n])
                new_lines += [line_]
                n=n+1

    data = ''.join(new_lines)
    with open(grofile, 'w') as f:    
        f.write(data)





    