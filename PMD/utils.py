import numpy as np
import mdtraj as md
import os

def get_atom_ids(pdb_name,res_ids):
    '''
    Takes pdb file name, and res_ids as input and returns the atoms ids
    Input:
        pdb_name:   string
        res_ids:    list of integers
    Return:
        atom_ids:   list of integers
    '''
    pdb = md.load_pdb(pdb_name)
    top = pdb.topology
    atom_ids = list()
    for res_id in res_ids:
        atom_ids.extend(top.select('resid '+str(res_id)))
    return atom_ids


def write_plumed_file(atom_ids,frequency,force,index,direction):
    '''
    write plumed input file, use BIAS_VALUE as BIAS potential
    Input:
        atom_ids:   list of integers
        frequency:  float (vibrational frequency)
        force:      float (transformed later to pN)
        direction:  str 'x', 'y' or 'z'
    Create new plumed.in No Return
    '''

    file_name = 'pump_'+str(index)+'_'+direction+'.in'
    atoms = [str(atom_id) for atom_id in atom_ids]
    atoms_str = ','.join(atoms)
    k = str(round(force/1.661,2))
    omega = str(round(frequency*6.28,2))
    f = open(file_name,'w')
    f.write('t: TIME\n')
    f.write('c: COM ATOMS='+atoms_str+'\n')
    f.write('k: MATHEVAL ARG=t VAR=t FUNC='+k+'*'+'cos(t*'+omega+') PERIODIC=NO\n')
    f.write('P: POSITION ATOM=c\n')
    f.write('V: MATHEVAL ARG=P.'+direction+',k VAR=x,y FUNC=x*y\n')
    f.write('BV: BIASVALUE ARG=V')
    f.close()


def write_plumed_files(atom_ids,frequency,force,index):
    '''
    write plumed input files for 3 directions
    '''
    write_plumed_file(atom_ids,frequency,force,index,'x')
    write_plumed_file(atom_ids,frequency,force,index,'y')
    write_plumed_file(atom_ids,frequency,force,index,'z')

def run_simulation(plumed_in,GPU,direction,index,cuda='0'):
    '''
    run simulations for pumped MD
    '''
    if os.path.exists('prod.in') and os.path.exists(plumed_in):
        output = 'prod'+str(index)+'_'+direction
        os.sys('export CUDA_VISIBLE_DEVICES='+cuda)
        out_file = output+'.out'
        mdcrd_file = output+'.mdcrd'
        mdvel_file = output+'.mdvel'
        run_command = 'pmemd.cuda -O -i prod.in -o '+out_file+' -p mol.prmtop -c equil.rst -x '+mdcrd_file+' -v '+mdvel_file
        os.sys(run_command)
    else:
        print('No md input files and please check your folder')


