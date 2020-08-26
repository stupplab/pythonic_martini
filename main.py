"""Python functions to create and simulate PAs
using GROMACS-5
Maintainer: Mayank Agrawal
"""



import numpy as np
import os
import json



name ='PA' # generic name for a peptide amphiphile


def gen_PA():
    # Write a gen_PA.pgn file specific to PA_name

    PA_name = 'PA_specific'

    data = [
    'package require psfgen',
    'psfcontext reset',
    'topology top_all36_prot_forPAs.rtf',
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
    'writepdb PA_aa.pdb',
    'exit'
    ]

    f = open('gen_PA.pgn', 'wt')
    f.write('\n'.join(data))
    f.close



def generic_to_specific_PA(PA_seq):
    """Write up a random PA sequence : PA1.pdb
    Exchange A,B,C... in PA_generic.pdb to valid residues
    """
    
    num_res = len(PA_seq)

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
               'Y','Z']
    
    with open('PA_generic.pdb', 'r') as fin:
        data = ''
        for line in fin:
            if 3*letters[num_res] in line:
                data = data + fin.readlines()[-1]
                break
            for i in range(num_res):
                if 3*letters[i] in line:
                    line = line.replace(3*letters[i], res_dict[PA_seq[i]])
                    break
            data = data + line

    with open('PA_specific.pdb', 'w') as fout:
        fout.write(data)



def make_aa_pdb(pep_seq):
    # Make the all-atom pdb file for the pep_seq
    # using vmd
    generic_to_specific_PA(pep_seq)
    gen_PA()
    os.system('vmd -dispdev text -e gen_PA.pgn')





def create_CGfiles_using_martinizepy(Ctermini_type, set_charge):
    """ Run martinize.py to create .top, .pdb, .itp files. 
    Ctermini_type: NH2 vs OH decide amide vs carboxylate (charged) termini
    Simplistic modification of .itp to put <charge> (is integer) 
    amount of charge on the farthest-from-hydrophobic side of PA.
    Raises error if putting charge is not plausible.
    NOTE: requires the heading [ atoms ] to be written as such 
    and rest of the lines to be contiguous - no linebreak
    """

    os.system('python2 martinize.py -f %s_aa.pdb \
        -o %s.top -x %s.pdb -name %s -ff martini22 \
        -nt \
        -ss CCCCCCCCCCCC '%(name,name,name,name))


    # Collect lines defining atoms
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
        
    

    # Modify lines_atoms as per Ctermini
    charged_thusfar = 0
    if Ctermini_type.upper() == 'OH':
        for i in range(len(lines_atoms))[::-1]:
            if 'BB' in lines_atoms[i]:
                lines_atoms[i]   = lines_atoms[i].replace(' 0.0', '-1.0')
                lines_atoms[i]   = lines_atoms[i].replace('P5', 'Qa')    
                charged_thusfar += -1
                break


    # modify charge of side chains,
    # CURRENTLY only neutralizes if Qd SC is found (deprotonation)
    neutralize_ahead = False
    if set_charge < 0: # deprotonation
        for i in range(len(lines_atoms))[::-1]:
            if charged_thusfar == set_charge:
                neutralize_ahead = True
            
            if ('SC' in lines_atoms[i]) and ('-1.0' in lines_atoms[i]):
                if 'Qa' not in lines_atoms[i]:
                    raise RuntimeError('-1.0 charge without Qa bead is found')
                if neutralize_ahead:
                    lines_atoms[i]   = lines_atoms[i].replace('-1.0', ' 0.0')
                    lines_atoms[i]   = lines_atoms[i].replace('Qa', 'P1')
                else:
                    charged_thusfar += -1

            if ('SC' in lines_atoms[i]) and (' 1.0' in lines_atoms[i]):
                if 'Qd' not in lines_atoms[i]:
                    raise RuntimeError('1.0 charge without Qd bead is found')
                lines_atoms[i]   = lines_atoms[i].replace('1.0', ' 0.0')
                lines_atoms[i]   = lines_atoms[i].replace('Qd', 'P1')

        if charged_thusfar != set_charge:
            raise ValueError('Peptide sequence could not be used to achieve set_charge')

    elif set_charge == 0: # protonation-deprotonation
        if Ctermini_type == 'OH':
            raise ValueError('Protonation after deprotonation does not make sense')
        
        for i in range(len(lines_atoms))[::-1]:
            if ('SC' in lines_atoms[i]) and ('-1.0' in lines_atoms[i]):
                if 'Qa' not in lines_atoms[i]:
                    raise RuntimeError('-1.0 charge without Qa bead is found')
                lines_atoms[i]   = lines_atoms[i].replace('-1.0', ' 0.0')
                lines_atoms[i]   = lines_atoms[i].replace('Qa', 'P1')

            if ('SC' in lines_atoms[i]) and (' 1.0' in lines_atoms[i]):
                if 'Qd' not in lines_atoms[i]:
                    raise RuntimeError('1.0 charge without Qd bead is found')
                lines_atoms[i]   = lines_atoms[i].replace('1.0', ' 0.0')
                lines_atoms[i]   = lines_atoms[i].replace('Qd', 'P1')
    
    elif set_charge > 0: # protonation
        if Ctermini_type == 'OH':
            raise ValueError('Protonation after deprotonation does not make sense')

        for i in range(len(lines_atoms))[::-1]:
            if charged_thusfar == set_charge:
                neutralize_ahead = True

            if ('SC' in lines_atoms[i]) and ('-1.0' in lines_atoms[i]):
                if 'Qa' not in lines_atoms[i]:
                    raise RuntimeError('-1.0 charge without Qa bead is found')
                lines_atoms[i]   = lines_atoms[i].replace('-1.0', ' 0.0')
                lines_atoms[i]   = lines_atoms[i].replace('Qa', 'P1')
                
            if ('SC' in lines_atoms[i]) and (' 1.0' in lines_atoms[i]):
                if 'Qd' not in lines_atoms[i]:
                    raise RuntimeError('1.0 charge without Qd bead is found')
                if neutralize_ahead:
                    lines_atoms[i]   = lines_atoms[i].replace('1.0', ' 0.0')
                    lines_atoms[i]   = lines_atoms[i].replace('Qd', 'P1')
                else:
                    charged_thusfar += 1
        
        if charged_thusfar != set_charge:
            raise ValueError('Peptide sequence could not be used to achieve set_charge')


    data_new = ''
    for line in data[:break1]:
        data_new += line
    for line in lines_atoms:
        data_new += line
    for line in data[break2:]:
        data_new += line
    
    
    with open('%s.itp'%name, 'w') as f:
        f.write(data_new)





def create_simulation_box(Lx,Ly,Lz,nmol):
    """Create simulation box
    Lx,Ly,Lz are box dimensions
    nmol are number of molecules inserted in the box,
    .pdb file created by martinizepy contains one such molecule
    """

    os.system('gmx insert-molecules \
        -box %s %s %s \
        -nmol %s \
        -ci %s.pdb \
        -radius 0.4 \
        -o %s_box.gro &> out.log'%(Lx,Ly,Lz,nmol,name,name))

    
    # Update number of PA molecules in .top file
    with open('out.log', 'r') as f:
        for line in f:
            if 'Added' in line:
                nmol = int(line.split()[1])
    with open('%s.top'%name, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]=='%s'%name:
                data+='%s       %s'%(name,nmol)
            else:
                data+=line
    with open('%s.top'%name, 'w') as f:
        f.write(data)




def solvate():
    water = 'water-12.5'
    os.system('gmx solvate \
        -cp %s_box.gro \
        -cs %s.gro \
        -radius 0.21 -o %s_water.gro &> out.log'%(name,water,name))

    # Append/Update number of water molecules in .top file
    with open('out.log', 'r') as f:
        for line in f:
            if 'Number of solvent molecules:' in line:
                nmol_W = int(line.split()[-1])
    added=False
    with open('%s.top'%name, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]=='W':
                data+='W        %s'%(nmol_W)
                added=True
            else:
                data+=line
    if not added:
        with open('%s.top'%name, 'a') as f:
            f.write('\nW        %s'%nmol_W)




def add_ions():
    # neutralize the PA_water.gro and update .top
    # TODO: Extra ions may be added to increase ion conc.

    """Old script that replaces water with ions 
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
    with open('%s.itp'%name, 'r') as f:
        charge = 0
        for line in f:
            if ('SC' in line) or ('BB' in line):
                if '-1.0' in line:
                    charge += -1
                elif '1.0' in line:
                    charge += 1
    
    # read the number of molecules added
    with open('%s.top'%name, 'r') as f:
        for line in f:
            if line.split() != []:
                if line.split()[0]=='PA':
                    nmol = int(line.split()[1])
                    

    n_ions = int(np.abs(nmol*charge))
    
    if charge < 0:
        # Add NA ions to neutralize the solution
        os.system('gmx insert-molecules \
            -f %s_water.gro \
            -nmol %s \
            -ci ion_NA.pdb \
            -radius 0.21 \
            -o %s_water.gro'%(name,n_ions,name))

    elif charge > 0:
        # Add CL ions to neutralize the solution
        os.system('gmx insert-molecules \
            -f %s_water.gro \
            -nmol %s \
            -ci ion_CL.pdb \
            -radius 0.21 \
            -o %s_water.gro'%(name,n_ions,name))

    if charge!=0:
        # include "martini_v2.0_ions.itp" in PA.top
        with open('%s.top'%name, 'r') as f:
            data = ''
            for line in f:
                if '#include "martini_v2.2.itp"' in line:
                    data+=line+'#include "martini_v2.0_ions.itp" \n'
                else:
                    data+=line
        with open('%s.top'%name, 'w') as f:
            f.write(data)


        # Append/Update number of ions in .top file
        if charge < 0:
            ion='NA+'
        else:
            ion='CL-'
        added=False
        with open('%s.top'%name, 'r') as f:
            data = ''
            for line in f:
                if len(line.split())!=0 and line.split()[0]==ion:
                    data+='%s      %s'%(ion,n_ions)
                    added=True
                else:
                    data+=line
        if not added:
            with open('%s.top'%name, 'a') as f:
                f.write('\n%s      %s'%(ion,n_ions))



def minimization():
    # Minimization
    os.system('gmx grompp \
        -f PA_water_min.mdp \
        -p %s.top \
        -c %s_water.gro \
        -o %s_water_min.tpr -maxwarn 1'%(name,name,name))
    os.system('gmx mdrun \
        -deffnm %s_water_min -v'%name)

    

def prepare_eq_tpr():
    # create PA_eq.tpr for mdrun
    os.system('gmx grompp \
            -f PA_water_eq.mdp \
            -p %s.top \
            -c %s_water_min.gro \
            -o %s_water_eq.tpr -maxwarn 2'%(name,name,name))



def equilibration(mpi=False):
    ## Equilibration
    ## Also prepare tpr
    if not mpi:
        os.system('gmx mdrun \
            -deffnm %s_water_eq -v -cpt 1'%name)
    else:
        os.system('mpirun -np 1 gmx_mpi mdrun \
            -deffnm %s_water_eq -v'%name)





def main(pep_seq, Lx,Ly,Lz, nmol, charge, Ctermini_type):
    
    make_aa_pdb(pep_seq)
    create_CGfiles_using_martinizepy(Ctermini_type, charge)
    create_simulation_box(Lx,Ly,Lz,nmol)
    solvate()
    add_ions()
    minimization()
    
    # prepare_eq_tpr()
    # equilibration(mpi=False)





if __name__=='__main__':
    
    if not os.path.isfile('signac_statepoint.json'):
        raise FileNotFoundError('signac_statepoint.json file needs to be created with all the parameters required here.\
            Example file:            \n \
            {                        \n \
              "pep_seq": "VEVE",     \n \
              "nmol": 100,           \n \
              "Lx": 12.5,            \n \
              "Ly": 12.5,            \n \
              "Lz": 12.5,            \n \
              "Ctermini_type": "OH", \n \
              "charge": -1           \n \
            }')

    
    with open('signac_statepoint.json', 'r') as f:
        sp = json.load(f)
        
        pep_seq       = sp['pep_seq']
        nmol          = sp['nmol']
        Lx            = sp['Lx']
        Ly            = sp['Ly']
        Lz            = sp['Lz']
        Ctermini_type = sp['Ctermini_type']
        charge        = sp['charge']
        
        main(pep_seq, Lx,Ly,Lz, nmol, charge, Ctermini_type)

    





