import numpy as np
import mdtraj as md
from np.fft import rfft, rfftfreq
import utils
import os

def calc_com(index,direction):
    '''
       load traj without water and ions
       and return the com traj of each residue
    '''
    traj_file = 'prod'+str(index)
    if index>1:
        traj_file = traj_file+'_'+direction
    traj_file + traj_file+'_nw.xtc'
    top_file = 'mol_nw.prmtop'
    traj = md.load(traj_file,top=top_file)
    top = traj.topology
    N = traj.n_residues
    com_data = list()
    for i in range(N):
        atom_ids = top.select('resid '+str(i))
        temp_traj = traj.atom_slice(atom_ids)
        com_data.append(md.compute_center_of_mass(temp_traj))
        del temp_traj
    return com_data


def calc_spectrum(frequency,range=0,index,direction):
    '''
    Claculate the baseline of unbiased simulation
    Inputs:
        frequency:  float the frequency of oscilation force
        range:      int take +-range number of frequencies in FFT for consideration default to be 0
    Output:
        return the normalized FFT as the function of residue index for the given frequency plus a small range
    '''


def calc_fft_intensity(data,frequency,range):
    '''
    Give a sequene of data points, calculate the FFT intensity of a particular frequency
    Input:
        data,       numpy array with length N
        frequency,  float frequency that calculate intensity at
        range,      int +- number of frequencies used to calculate intensity
    Returns:
        the normalized intensity at given frequency
    
    N = len(data)
    yf = rfft(data)
    xf = rfftfreq() 
    '''
    return None

def write_cpptraj_vac_input_file(pdb_file,crd_file,vel_file,top_file):
    '''
    write cpptraj input file to calculate VAC for each residue
    Input:
        crd_file:   crd file name
        vel_file;   vel file name
        top_file:   topology file name
    '''
    res_num = utils.get_res_num(pdb_file)
    file_name='ctj_vac.in'
    f = open(file_name,'w')
    f.write('parm '+top_file+'\n')
    f.write('trajin '+crd_file+' mdvel '+vel_file+'\n')
    for i in range(res_num):
        atom_ids = utils.get_atom_ids(pdb_file,i+1)
        atom_ids_str_list = [str(id+1) for id in atom_ids]
        atom_mask = ','.join(atom_ids_str_list)
        out_file = 'res_'+str(i+1)+'_vac.out'
        f.write('velocityautocorr @'+atom_mask+' out '+out_file+' tstep 0.01 norm\n')
    f.close()
    return None
