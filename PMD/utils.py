import numpy as np
import mdtraj as md
import os

def get_atom_ids(pdb_name,res_ids):
    '''
    Takes pdb file name, and res_ids as input and returns the atoms ids
    Request: prepare the pdb that has added H atoms
    Inputs:
        pdb_name:   string
        res_ids:    list of integers
    Return:
        atom_ids:   list of integers
    '''
    pdb = md.load_pdb(pdb_name)
    top = pdb.topology
    atom_ids = list()
    if isinstance(res_ids,int):
        atom_ids=top.select('resid '+str(res_ids))
    else:
        for res_id in res_ids:
            atom_ids.extend(top.select('resid '+str(res_id)))
    return atom_ids


def write_plumed_file(atom_ids,frequency,force,index,direction):
    '''
    write plumed input file, use BIAS_VALUE as BIAS potential
    Inputs:
        atom_ids:   list of integers
        frequency:  float (vibrational frequency)
        force:      float (transformed later to pN)
        direction:  str 'x', 'y' or 'z'
    Create new plumed.in No Return
    '''

    file_name = 'pump'+str(index)+'_'+direction+'.dat'
    atoms = [str(atom_id) for atom_id in atom_ids]
    atoms_str = ','.join(atoms)
    k = str(round(force/1.661,2))
    omega = str(round(frequency*6.28,2))
    f = open(file_name,'w')
    f.write('t: TIME\n')
    f.write('c: COM ATOMS='+atoms_str+'\n')
    f.write('k: MATHEVAL ARG=t VAR=t FUNC='+k+'*'+'cos(t*'+omega+') PERIODIC=NO\n')
    f.write('P: POSITION ATOM=c\n')
    f.write('V: MATHEVAL ARG=P.'+direction+',k VAR=x,y FUNC=x*y PERIODIC=NO\n')
    f.write('BV: BIASVALUE ARG=V')
    f.close()


def write_plumed_files(atom_ids,frequency,force,index):
    '''
    write plumed input files for 3 directions
    Inputs:
        atom_ids:   list of integers
        frequency:  float (see write_plumed_file)
        force:      float (see write_plumed_file)
        direction:  str   (see write_plumed_file)
    '''
    write_plumed_file(atom_ids,frequency,force,index,'x')
    write_plumed_file(atom_ids,frequency,force,index,'y')
    write_plumed_file(atom_ids,frequency,force,index,'z')

def run_simulation(direction,index,time,cuda='0'):
    '''
    run simulations for pumped MD
    Inputs:
        direction:  str
        index:      int
        cuda:       str (optional define the GPU name to use)
        time:       int in ps
    '''
    input_file = write_prod_input_file(plumed_in,index,direction,time)
    output = 'prod'+str(index)
    if index>0:
        output = output+'_'+direction
    cuda_command = 'export CUDA_VISIBLE_DEVICES='+cuda
    os.system(cuda_command)
    out_file = output+'.out'
    mdcrd_file = output+'.mdcrd'
    #mdvel_file = output+'.mdvel'
    run_command = 'pmemd.cuda -O -i '+input_file+' -o '+out_file+' -p mol.prmtop -c equil.rst -x '+mdcrd_file
    os.system(run_command)

def write_prod_input_file(direction,index,time):
    '''
    write prod input files
    Inputs:
        time:   int in ps
        index:  if index=0 it is unbiased simulation
    return input file name
    '''
    nsteps = str(int(1000*time))
    file_name = 'prod'+str(index)
    if index>0:
        file_name = file_name+'_'+direction
    file_name = file_name+'.in'
    f = open(file_name,'w')
    f.write(file_name+'\n')
    f.write('&cntrl\n')
    f.write('  imin=0,irest=1,ntx=5,\n  nstlim='+nsteps+',dt=0.001,\n  ntc=2,ntf=2,\n')
    f.write('  cut=8.0, ntb=2, ntp=1, taup=2.0,\n')
    f.write('  ntpr=5000,ntwx=10,\n')
    f.write('  ntt=3,gamma_ln=2.0,temp0=300.0,ig=-1\n')
    if index>0:
        plumed_in = 'pump'+str(index)+'_'+direction+'.dat'
        f.write(',\n  plumed=1,plumedfile=\''+plumed_in+'\'\n/\n')
    else:
        f.write('/\n')
    f.close()
    return file_name

def write_cpptraj_input_file(index,direction):
    '''
    write cpptraj inpute file for extract waters and Na+ and CL-
    return the file name
    '''
    file_name='ctj'+str(index)
    trajin = 'trajin prod'+str(index)
    trajout = 'trajout prod'+str(index)
    if index>0:
        file_name = file_name+'_'+direction
        trajin = trajin+'_'+direction
        trajout = trajout+'_'+direction
    file_name=file_name+'.in'
    trajin = trajin+'.mdcrd'
    trajout = trajout+'_nw.xtc'
    f = open(file_name,'w')
    f.write('parm mol.prmtop\n')
    f.write(trajin+'\n')
    f.write('strip :WAT\n')
    f.write('strip :Na+\n')
    f.write('strip :Cl-\n')
    f.write(trajout+'\n')
    f.close()
    return file_name 

def run_cpptraj(index,direction):
    '''
    run cpptraj to remove waters and ions of a trajectory and save the file
    Inputs:
        ctj_in: str name of cpptraj input file
    '''
    ctj_in = write_cpptraj_input_file(index,direction)
    os.system('cpptraj -i '+ctj_in)
    return None

def strip_topology_wat():
    '''
    strip water and ions from topology file
    '''
    if os.path.exists('mol.prmtop'):
        f = open('ctj_top.in','w')
        f.write('parm mol.prmtop\n')
        f.write('parmstrip :WAT\n')
        f.write('parmstrip :Na+\n')
        f.write('parmstrip :Cl-\n')
        f.write('parmwrite out mol_nw.prmtop\n')
        f.close()
        os.system('cpptraj -i ctj_top.in')
    else:
        print('No mol.prmtop file, maybe in the wrong directory.')
    return None
