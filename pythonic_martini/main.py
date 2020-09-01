"""Python functions to create and simulate PAs
tested using GROMACS-5
Maintainer: Mayank Agrawal
"""



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
    """Write up a random PA sequence : <name>.pdb
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
               'Y','Z'] # generic residue name in PA_generic.pdb
    
    with open('%s/PA_generic.pdb'%this_path, 'r') as fin:
        data = ''
        for line in fin:
            if 3*letters[num_res] in line:
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








def make_aa_pdb(PA_seq, name):
    """Make the all-atom pdb file <name>_aa.pdb for the PA_seq using vmd
    PA_seq is in string format. Example: C12VVAAEE. 
    Digits are used only for the alkyl chain.
    """
    generic_to_specific_PA(PA_seq.upper(), name)
    gen_PA(name)
    os.system('vmd -dispdev text -e gen_%s.pgn'%name)





def create_CGfiles_using_martinizepy(Ctermini_type, res_charge=[], name='pep'):
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

    os.system('python2 %s/martinize.py -f %s_aa.pdb \
        -o %s.top -x %s.pdb -name %s -ff martini22 \
        -nt \
        -ss CCCCCCCCCCCC '%(this_path,name,name,name,name))


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

    os.system('gmx insert-molecules \
        -box %s %s %s \
        -nmol %s \
        -ci %s.pdb \
        -radius %s \
        -o %s &> out.log'%(Lx,Ly,Lz,nmol,molname,vdwradius,boxfilename))

    
    # Actual nmol added
    with open('out.log', 'r') as f:
        for line in f:
            if 'Added' in line:
                nmol = int(line.split()[1])
    
    # read <molname>.top
    with open('%s.top'%molname, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]=='%s'%molname:
                data+='%s'%molname + ' '*(9-len(molname)) + '%s\n'%nmol
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
    os.system('gmx insert-molecules \
            -f %s \
            -nmol %s \
            -ci %s.pdb \
            -radius %s \
            -o %s &> out.log'%(inwhichfilename,nmol,molname,vdwradius,outfilename))

    # Actual nmol added
    with open('out.log', 'r') as f:
        for line in f:
            if 'Added' in line:
                nmol = int(line.split()[1])


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
                data+=molname + ' '*(9-len(molname)) + '%s\n'%nmol
                mol_added = True
            else:
                data+=line

        if not mol_added:
            data += molname + ' '*(9-len(molname)) + '%s\n'%nmol
    
    with open('%s'%topfilename, 'w') as f:
        f.write(data)




def insert_water(inwhichfilename, volume, vdwradius, outfilename, topfilename):
    
    # One martini water bead is 72 g/mol (4 times the H2O)
    num_water = int( 1000 *volume / (1.66*72) )
    molname = 'W'

    os.system('gmx insert-molecules \
            -f %s \
            -nmol %s \
            -ci %s/W.pdb \
            -radius %s \
            -o %s &> out.log'%(inwhichfilename,num_water,this_path,vdwradius,outfilename))

    # Actual nmol added
    with open('out.log', 'r') as f:
        for line in f:
            if 'Added' in line:
                nmol = int(line.split()[1])


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
    os.system('gmx solvate \
        -cp %s \
        -cs %s/%s.gro \
        -radius %s -o %s &> out.log'%(inwhichfilename,this_path,water,vdwradius,outfilename))

    
    # Append/Update number of water molecules in .top file
    with open('out.log', 'r') as f:
        for line in f:
            if 'Number of solvent molecules:' in line:
                nmol_W = int(line.split()[-1])
    
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




def neutralize_system(inwhichfilename, molname, vdwradius, outfilename, topfilename):
    # neutralize the <molname> using the charge from <molname>.itp and update <topfilename>

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
        os.system('gmx insert-molecules \
            -f %s \
            -nmol %s \
            -ci %s/ion_NA.pdb \
            -radius %s \
            -o %s'%(inwhichfilename,nions,this_path,vdwradius,outfilename))

    elif charge > 0:
        # Add CL ions to neutralize the solution
        os.system('gmx insert-molecules \
            -f %s \
            -nmol %s \
            -ci %s/ion_CL.pdb \
            -radius %s \
            -o %s'%(inwhichfilename,nions,this_path,vdwradius,outfilename))

    if charge!=0:
        os.system('cp %s/%s ./'%(this_path,martini_ionitp))
        
        # include martini_ionitp in topfilename if not already included
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
            


        # Append/Update number of ions in .top file
        # adds the ion amount already present
        if charge < 0:
            ion='NA+'
        else:
            ion='CL-'
        added=False
        with open(topfilename, 'r') as f:
            data = ''
            for line in f:
                if len(line.split())!=0 and line.split()[0]==ion:
                    data+=ion+ ' '*(9-len(ion)) + '%s\n'%(nions+int(line.split()[1]))
                    added=True
                else:
                    data+=line
        if not added:
            data += '%s'%ion + ' '*(9-len(ion)) + '%s\n'%nions
        with open(topfilename, 'w') as f:
            f.write(data)




def add_ions(inwhichfilename, nions, vdwradius, outfilename, topfilename):
    # <nions> number of NA and CL ions are added to increase ionic concentration    
    
    if nions==0:
        warnings.warm('Zero ions are asked to add. Doing nothing.')
        return
        
    os.system('cp %s/%s ./'%(this_path,martini_ionitp))

    # add extra ions
    os.system('gmx insert-molecules \
        -f %s \
        -nmol %s \
        -ci %s/ion_NA.pdb \
        -radius %s \
        -o %s'%(inwhichfilename,nions,this_path,vdwradius,outfilename))

    os.system('gmx insert-molecules \
        -f %s \
        -nmol %s \
        -ci %s/ion_CL.pdb \
        -radius %s \
        -o %s'%(inwhichfilename,nions,this_path,vdwradius,outfilename))


    ## update topfilename

    # include martini_ionitp in topfilename if not already included
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

    # Append/Update number of ions in .top file
    # adds the ion amount already present
    
    ion='NA+'
    added=False
    with open(topfilename, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]==ion:
                data+=ion+ ' '*(9-len(ion)) + '%s\n'%(nions+int(line.split()[1]))
                added=True
            else:
                data+=line
    if not added:
        data += '%s'%ion + ' '*(9-len(ion)) + '%s\n'%nions
    with open(topfilename, 'w') as f:
        f.write(data)

    ion='CL-'
    added=False
    with open(topfilename, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]==ion:
                data+=ion+ ' '*(9-len(ion)) + '%s\n'%(nions+int(line.split()[1]))
                added=True
            else:
                data+=line
    if not added:
        data += '%s'%ion + ' '*(9-len(ion)) + '%s\n'%nions
    with open(topfilename, 'w') as f:
        f.write(data)






def minimization(mdpfilename, topfilename, grofilename, tprfilename):
    # Minimization
    os.system('gmx grompp \
        -f %s \
        -p %s \
        -c %s \
        -o %s -maxwarn 1'%(mdpfilename,topfilename,grofilename,tprfilename))
    
    tprname = os.path.basename(tprfilename).strip('.tpr')
    os.system('gmx mdrun \
        -deffnm %s -v'%tprname)

    

def prepare_eq_tpr(mdpfilename, topfilename, grofilename, tprfilename):
    # create .tpr for mdrun
    
    os.system('gmx grompp \
        -f %s \
        -p %s \
        -c %s \
        -o %s -maxwarn 1'%(mdpfilename,topfilename,grofilename,tprfilename))



def equilibration(tprfilename, mpi=False):
    ## Equilibration

    tprname = os.path.basename(tprfilename).strip('.tpr')
    ## Also prepare tpr
    if not mpi:
        os.system('gmx mdrun \
            -deffnm %s -v -cpt 30'%tprname)
    else:
        os.system('mpirun -np 1 gmx_mpi mdrun \
            -deffnm %s -v'%tprname)








if __name__=='__main__':
    
    print('main.py should not be run directly, it to be used in a script. See example_script.py')

    





