"""Python functions to create and simulate PAs
tested using GROMACS-5
Maintainer: Mayank Agrawal



Additional info;

out.log file, you find here, is only used for dumping the output of some commands
"""


import subprocess
import numpy as np
import os
import json
import warnings


this_path = os.path.dirname(os.path.abspath(__file__))
martini_itp    = "martini_v2.2.itp"
martini_ionitp = "martini_v2.0_ions.itp"


amino_acid_1to3_lettercode = {'A':"ALA",'C':"CYS",'D':"ASP",'E':"GLU",'F':"PHE", 
                'G':"GLY",'H':"HIS",'I':"ILE",'K':"LYS",'L':"LEU", 
                'M':"MET",'N':"ASN",'P':"PRO",'Q':"GLN",'R':"ARG", 
                'S':"SER",'T':"THR",'V':"VAL",'W':"TRP",'Y':"TYR"}



def gen_PA(name):
    # Write a gen_PA.pgn file specific to PA_name

    PA_name = '%s_specific'%name

    data = [
    'package require psfgen',
    'psfcontext reset',
    'topology %s/top_all36_prot_forPAs.rtf'%this_path,
    'pdbalias residue HIS HSD',
    'segment P001 {',
    'pdb %s.pdb'%PA_name,
    'first none',
    '#last CT2',
    'auto angles dihedrals',
    '}',
    'coordpdb %s.pdb P001'%PA_name,
    '#last CTER',
    'regenerate angles dihedrals',
    'guesscoord',
    'writepdb %s_aa.pdb'%name,
    'exit'
    ]

    f = open('gen_%s.pgn'%name, 'wt')
    f.write('\n'.join(data))
    f.close



def generic_to_specific_PA(PA_seq, name):
    """Write up a specific PA sequence PA_seq : <name>_specific.pdb
    Use the amount of alkyl chain asked for. Currently the maximum is C16.
    Exchange A,B,C... in PA_generic.pdb to valid residues
    """
    
    if (PA_seq[0] == 'C') and PA_seq[1].isdigit():
        if PA_seq[2].isdigit():
            num_alkylC = int(PA_seq[1:3])
            pep_seq = PA_seq[3:]
        else:
            num_alkylC = int(PA_seq[1])
            pep_seq = PA_seq[2:]
    else: # no alkyl
        pep_seq = PA_seq
        num_alkylC = 0

    num_res = len(pep_seq)

    
    res_dict = {'A':"ALA",'C':"CYS",'D':"ASP",'E':"GLU",'F':"PHE", 
                'G':"GLY",'H':"HIS",'I':"ILE",'K':"LYS",'L':"LEU", 
                'M':"MET",'N':"ASN",'P':"PRO",'Q':"GLN",'R':"ARG", 
                'S':"SER",'T':"THR",'V':"VAL",'W':"TRP",'Y':"TYR"}

    # residues = ['ALA','DAL','ARG','ASN','ASP','CYS',
    #             'GLN','GLU','DGL','GLY','HSD','HSE',
    #             'HSP','DHD','DHE','DHP','ILE','LEU',
    #             'LYS','DLY','MET','PHE','PRO','SER',
    #             'DSE','THR','TRP','TYR','VAL','DVA',
    #             'PAM','ALAD']

    letters = ['A','B','C','D','E','F','G','H','I','J','K','L',
               'M','N','O','P','Q','R','S','T','U','V','W','X',
               'Y','Z'] + [
              'A','B','C','D','E','F','G','H','I','J','K','L',
               'M','N','O','P','Q','R','S','T','U','V','W','X',
               'Y','Z']
               # generic residue name in PA_generic.pdb
    
    with open('%s/PA_generic.pdb'%this_path, 'r') as fin:
        data = ''
        for line in fin:
            if (3*letters[num_res] in line) or (letters[num_res]+'ZZ' in line):
                data = data +'END'
                break
            if 'PAM' in line:
                if int(line.split()[2][1:]) > num_alkylC:
                    continue
            for i in range(num_res):
                if 3*letters[i] in line:
                    line = line.replace(3*letters[i], res_dict[pep_seq[i]])
                    break
            data = data + line

    with open('%s_specific.pdb'%name, 'w') as fout:
        fout.write(data)





def _actual_atoms_added(filename, keyword):
    """Count the number of lines with <keyword> in the <filename> - usually gro or pdb file
    """
    num_keyword = 0
    with open(filename, 'r') as f:
        for line in f:
            if keyword in line:
                num_keyword += 1

    return num_keyword




def _actual_molecules_added(filename, itpfilename, start_linenumber=0):
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





def make_aa_pdb(PA_seq, name):
    """Make the all-atom pdb file <name>_aa.pdb for the PA_seq using vmd
    PA_seq is in string format. Example: C12VVAAEE. 
    Digits are used only for the alkyl chain.
    """
    generic_to_specific_PA(PA_seq.upper(), name)
    gen_PA(name)
    os.system('vmd -dispdev text -e gen_%s.pgn'%name)

    # replace HSD back to HIS to be recognized by martinize.py
    with open('%s_aa.pdb'%name, 'r') as f:
        data = f.read().replace(' HSD ', ' HIS ')

    with open('%s_aa.pdb'%name, 'w') as f:
        f.write(data)




def create_CGfiles_using_martinizepy(Ctermini_type, res_charge=[], name='pep', ss='C'*20):
    """ Run martinize.py to create <name>.top, <name>.pdb, <name>.itp files. 
    Ctermini_type: NH2 vs OH decide amide vs carboxylate (charged) termini
    Modifies .itp to set <res_charge> as asked.
    
    res_charge is the list of (residue,charge) tuples 
    in order to how they appear in the <name>.itp and the sequence given to make_aa_pdb.
    Only K,L,D,E residues are looked at for charging or discharging, rest are ignored.
    Extra residues (e.g. 2 Ks are defined for the peptide with 1 K), if defined, are also ignored
    
    Raises error if putting charge is not plausible.
    NOTE: requires the heading [ atoms ] in .itp to be written as such 
    and rest of the lines to be contiguous - no linebreak
    """

    os.system('cp %s/%s ./'%(this_path,martini_itp))
    
    process = subprocess.run(f'python2 {this_path}/martinize.py -f {name}_aa.pdb \
        -o {name}.top -x {name}.pdb -name {name} -ff martini22 \
        -nt \
        -ss {ss} ', shell=True)
    process.check_returncode()
    
    
    # Collect lines defining atoms, breaks define start and end
    lines_atoms = []
    break1,break2 = None,None
    with open('%s.itp'%name, 'r') as f:
        data = f.readlines()
        start = False
        for i,line in enumerate(data):
            if '[ atoms ]' in line:
                start = True
                break1 = i+1
                continue
            if start:
                if line.split()==[]:
                    start = False
                    break2 = i
                    break
                lines_atoms = lines_atoms + [line]
        
    

    
    # The above code puts neutral P5 at C-terminus, which corresponds to the NH2 termini
    # Below code sets the Ctermini OH, at the last BB bead found
    charged_thusfar = 0
    if Ctermini_type.upper() == 'OH':
        for i in range(len(lines_atoms))[::-1]:
            if 'BB' in lines_atoms[i]:
                lines_atoms[i]   = lines_atoms[i].replace(' 0.0', '-1.0')
                lines_atoms[i]   = lines_atoms[i].replace('P5', 'Qa')    
                charged_thusfar += -1
                break


    # set <res_charge> in <name>.itp
    # Charged amino acids are: (+ive) K, R ; (-ive) D, E
    res_charge_SC = []
    for res,charge in res_charge:
        if res.upper() in 'KR':
            res_charge_SC += [[amino_acid_1to3_lettercode[res],charge,'SC2']]
        elif res.upper() in 'DE':
            res_charge_SC += [[amino_acid_1to3_lettercode[res],charge,'SC1']]

    k=0
    for i,line in enumerate(lines_atoms):
        if k==len(res_charge_SC):
            break
        res,charge,SC = res_charge_SC[k]
        if (res in line) and (SC in line):
            if (res == 'LYS'):
                if charge == 0:
                    lines_atoms[i]   = lines_atoms[i].replace(' 1.0', ' 0.0')
                    lines_atoms[i]   = lines_atoms[i].replace('Qd', 'P1')
                elif charge == 1:
                    lines_atoms[i]   = lines_atoms[i].replace(' 0.0', ' 1.0')
                    lines_atoms[i]   = lines_atoms[i].replace(line.split()[1], 'Qd')
                else:
                    raise ValueError('Charge other than 0 or 1 is set for LYS')

            elif (res == 'ARG'):
                if charge == 0:
                    lines_atoms[i]   = lines_atoms[i].replace(' 1.0', ' 0.0')
                    lines_atoms[i]   = lines_atoms[i].replace('Qd', 'P4')
                elif charge == 1:
                    lines_atoms[i]   = lines_atoms[i].replace(' 0.0', ' 1.0')
                    lines_atoms[i]   = lines_atoms[i].replace(line.split()[1], 'Qd')
                else:
                    raise ValueError('Charge other than 0 or 1 is set for ARG')
            
            elif (res == 'GLU'):
                if charge == 0:
                    lines_atoms[i]   = lines_atoms[i].replace('-1.0', ' 0.0')
                    lines_atoms[i]   = lines_atoms[i].replace('Qa', 'P1')
                elif charge == -1:
                    lines_atoms[i]   = lines_atoms[i].replace(' 0.0', '-1.0')
                    lines_atoms[i]   = lines_atoms[i].replace(line.split()[1], 'Qa')
                else:
                    raise ValueError('Charge other than 0 or -1 is set for GLU')
            
            elif (res == 'ASP'):
                if charge == 0:
                    lines_atoms[i]   = lines_atoms[i].replace('-1.0', ' 0.0')
                    lines_atoms[i]   = lines_atoms[i].replace('Qa', 'P3')
                elif charge == -1:
                    lines_atoms[i]   = lines_atoms[i].replace(' 0.0', '-1.0')
                    lines_atoms[i]   = lines_atoms[i].replace(line.split()[1], 'Qa')
                else:
                    raise ValueError('Charge other than 0 or 1 is set for LYS')

            k=k+1
            


    data_new = ''
    for line in data[:break1]:
        data_new += line
    for line in lines_atoms:
        data_new += line
    for line in data[break2:]:
        data_new += line
    

    with open('%s.itp'%name, 'w') as f:
        f.write(data_new)





def create_simulation_box(Lx,Ly,Lz,nmol,molname,vdwradius,boxfilename, topfilename=None):
    """Create simulation box
    Lx,Ly,Lz are box dimensions
    nmol are number of molecules inserted in the box,
    .pdb file created by martinizepy contains one such molecule
    """

    process = subprocess.run(f'gmx insert-molecules \
        -box {Lx} {Ly} {Lz} \
        -nmol {nmol} \
        -ci {molname}.pdb \
        -radius {vdwradius} \
        -o {boxfilename} &> out.log', shell=True)
    process.check_returncode()

    
    
    # Calculate actual nmol added
    actual_nmol = _actual_molecules_added(boxfilename, '%s.itp'%molname, start_linenumber=0)

    
    # read <molname>.top
    with open('%s.top'%molname, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]=='%s'%molname:
                data+='%s'%molname + ' '*(9-len(molname)) + '%s\n'%actual_nmol
            else:
                data+=line
        
    if topfilename==None:
        # Update number of PA molecules in <molname>.top file    
        with open('%s.top'%molname, 'w') as f:
            f.write(data)
    else:
        # Create <newtopfilename> file    
        with open('%s'%topfilename, 'w') as f:
            f.write(data)





def insert_molecules(inwhichfilename, molname, nmol, vdwradius, outfilename, topfilename, add_itp=True):
    """Insert molecules given in <what_filename> in the existing simulation box <inwhich_filename>
    Creates <outfilename> and update existing <topfilename>
    """
    
    with open(inwhichfilename, 'r') as f:
        num_lines = len(f.readlines())

    process = subprocess.run(f'gmx insert-molecules \
            -f {inwhichfilename} \
            -nmol {nmol} \
            -ci {molname}.pdb \
            -radius {vdwradius} \
            -o {outfilename} &> out.log', shell=True)
    process.check_returncode()

    
    # Calculate actual nmol added
    start_linenumber = num_lines-1
    actual_nmol = _actual_molecules_added(outfilename, '%s.itp'%molname, start_linenumber)


    # Update .top file to 
    # add <molname>.itp 
    # append/update number of <molname> molecules
    itp_added = not add_itp
    mol_added = False
    with open('%s'%topfilename, 'r') as f:
        data = ''
        for line in f:    
            if (not itp_added) and (len(line.split())==0): # add itp in the first available line
                data+='#include "%s.itp"\n'%molname +line
                itp_added = True
            if len(line.split())!=0 and line.split()[0]=='%s'%molname:
                data+=molname + ' '*(9-len(molname)) + '%s\n'%actual_nmol
                mol_added = True
            else:
                data+=line

        if not mol_added:
            data += molname + ' '*(9-len(molname)) + '%s\n'%actual_nmol
    
    with open('%s'%topfilename, 'w') as f:
        f.write(data)




def insert_water(inwhichfilename, volume, vdwradius, outfilename, topfilename):
    
    # One martini water bead is 72 g/mol (4 times the H2O)
    num_water = int( 1000 *volume / (1.66*72) )
    molname = 'W'

    process = subprocess.run(f'gmx insert-molecules \
            -f {inwhichfilename} \
            -nmol {num_water} \
            -ci {this_path}/W.pdb \
            -radius {vdwradius} \
            -o {outfilename} &> out.log', shell=True)
    process.check_returncode()

    
    # Actual W atoms added
    nmol = _actual_atoms_added(outfilename, 'W ')
    

    # Update .top file to 
    # append/update number of <molname> molecules
    mol_added = False
    with open('%s'%topfilename, 'r') as f:
        data = ''
        for line in f:    
            if len(line.split())!=0 and line.split()[0]=='%s'%molname:
                data+='%s'%molname + ' '*(9-len(molname)) + '%s\n'%nmol
                mol_added = True
            else:
                data+=line

        if not mol_added:
            data += '%s'%molname + ' '*(9-len(molname)) + '%s\n'%nmol
    
    with open('%s'%topfilename, 'w') as f:
        f.write(data)





def solvate(inwhichfilename, vdwradius, topfilename, outfilename):
    """Add water to the simulation box using a premade gro file - water-12.5.gro
    """

    water = 'water-12.5'
    process = subprocess.run(f'gmx solvate \
        -cp {inwhichfilename} \
        -cs {this_path}/{water}.gro \
        -radius {vdwradius} -o {outfilename} &> out.log')
    process.check_returncode()

    
    # Actual W atoms added
    nmol_W = _actual_atoms_added(outfilename, 'W ')
    

    added=False
    with open(topfilename, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]=='W':
                data+='W' + ' '*8 + '%s\n'%nmol_W
                added=True
            else:
                data+=line
        if not added:
            data+='W' + ' '*8 + '%s\n'%nmol_W
    

    with open(topfilename, 'w') as f:
        f.write(data)




def add_ions(inwhichfilename, ionname, nions, outfilename, topfilename):
    """
    Replaces <nions> number of W in <inwhichfilename> with ion <ionname>.
    Updates the topfile to correct number of W and ions

    ionname: Currently only accepts NA or CL

    This method i used by neutralize_system and increase_ionic_strength
    """

    if nions==0:
        warnings.warm('Zero %s ions are asked to add. Doing nothing.'%ionname)
        return

    
    if 'NA' in ionname.upper():
        ionname = 'NA'
        ion     = 'NA+'
    elif 'CL' in ionname.upper():
        ionname = 'CL'
        ion     = 'CL-'

    else:
        raise ValueError('%s not recognized'%ionname)

    

    # include martini_ionitp in topfilename if not already included
    os.system('cp %s/%s ./'%(this_path,martini_ionitp))

    with open(topfilename, 'r') as f:
        data = f.read()    
    
    if not '#include "%s"'%martini_ionitp in data:
        with open(topfilename, 'r') as f:
            data = ''
            for line in f:
                if '#include "%s"'%martini_itp in line:
                    data+=line+'#include "%s"\n'%martini_ionitp
                else:
                    data+=line
        with open(topfilename, 'w') as f:
            f.write(data)
    

    

    # replace last <nions> W atoms with <ionname>
    with open(inwhichfilename, 'r') as f:
        lines = f.readlines()
        
    added_ions = 0
    lines_reversed = lines[::-1]
    for i,line in enumerate(lines_reversed):
        if added_ions == nions:
            break
        if 'W ' in line:
            lines_reversed[i] = lines_reversed[i].replace(' '*(len(ionname)-1)+'W', ionname)
            lines_reversed[i] = lines_reversed[i].replace('W'+' '*(len(ion)-1), ion)
            added_ions += 1
    lines = lines_reversed[::-1]

    # write outfile
    with open(outfilename, 'w') as f:
        f.write(''.join(lines))



    # recalculate nions and W in the new file
    # recalculate W atoms added
    net_W = _actual_atoms_added(outfilename, 'W ')
    net_nions = _actual_atoms_added(outfilename, ion)



    # Append/Update number of ions in .top file
    added=False
    with open(topfilename, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]==ion:
                data+=ion+ ' '*(9-len(ion)) + '%s\n'%(net_nions)
                added=True
            elif len(line.split())!=0 and line.split()[0]=='W':
                data+='W'+ ' '*(9-len('W')) + '%s\n'%(net_W)
            else:
                data+=line
    if not added:
        data += '%s'%ion + ' '*(9-len(ion)) + '%s\n'%net_nions
    with open(topfilename, 'w') as f:
        f.write(data)





def neutralize_system(inwhichfilename, molname, outfilename, topfilename):
    """
    # neutralize the <molname> using the charge from <molname>.itp and update <topfilename>
    Replace W molecues in <inwhichfilename> and updates topfilename

    Old script that replaces water with ions 
    but requires command line input to choose the continuous solvent group.
    Hence, automation is difficult
    os.system('gmx grompp \
        -f peptide_water_min.mdp \
        -p %s.top \
        -c %s_water.gro \
        -o %s_genion.tpr'%(name,name,name))
    os.system('gmx genion \
        -s %s_genion.tpr \
        -pname NA+ -nname CL- -neutral \
        -p %s.top -o %s_water.gro && 3'%(name,name,name))
    """

    
    # read net charge on PA from.itp
    with open('%s.itp'%molname, 'r') as f:
        charge = 0
        for line in f:
            if ('SC' in line) or ('BB' in line):
                if '-1.0' in line:
                    charge += -1
                elif '1.0' in line:
                    charge += 1
    

    # read the number of molecules added
    with open(topfilename, 'r') as f:
        for line in f:
            if line.split() != []:
                if line.split()[0]==molname:
                    nmol = int(line.split()[1])
                    

    nions = int(np.abs(nmol*charge))
    

    if charge < 0:
        # Add NA ions to neutralize the solution
        add_ions(
            inwhichfilename,
            'NA',
            nions,
            outfilename,
            topfilename)


    elif charge > 0:
        # Add CL ions to neutralize the solution
        add_ions(
            inwhichfilename,
            'CL',
            nions,
            outfilename,
            topfilename)





def increase_ionic_strength(inwhichfilename, nions, outfilename, topfilename):
    # <nions> number of both NA and CL ions are added to increase ionic concentration
    
    if nions==0:
        warnings.warm('Zero ions are asked to add. Doing nothing.')
        return
        
    add_ions(
        inwhichfilename,
        'NA',
        nions,
        outfilename,
        topfilename)

    add_ions(
        inwhichfilename,
        'CL',
        nions,
        outfilename,
        topfilename)

    




def minimization(mdpfilename, topfilename, grofilename, tprfilename):
    # Minimization
    process = subprocess.run(f'gmx grompp \
        -f {mdpfilename} \
        -p {topfilename} \
        -c {grofilename} \
        -o {tprfilename} -maxwarn 1', shell=True)
    process.check_returncode()

    tprname = os.path.basename(tprfilename).replace('.tpr','')
    process = subprocess.run('gmx mdrun \
        -deffnm %s -v'%tprname, shell=True)
    process.check_returncode()    

    

def prepare_eq_tpr(mdpfilename, topfilename, grofilename, tprfilename):
    # create .tpr for mdrun
    
    process = subprocess.run(f'gmx grompp \
        -f {mdpfilename} \
        -p {topfilename} \
        -c {grofilename} \
        -o {tprfilename} -maxwarn 1', shell=True)
    process.check_returncode()    



def equilibration(tprfilename, mpi=False):
    ## Equilibration

    tprname = os.path.basename(tprfilename).replace('.tpr','')
    ## Also prepare tpr
    if not mpi:
        process = subprocess.run(f'gmx mdrun \
            -deffnm %s -v -cpt 30'%tprname, shell=True)
        process.check_returncode()

    else:
        process = subprocess.run(f'mpirun -np 1 gmx_mpi mdrun \
            -deffnm %s -v'%tprname, shell=True)
        process.check_returncode()







if __name__=='__main__':
    
    print('main.py should not be run directly, it to be used in a script. See example_script.py')

    





