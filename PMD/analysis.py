import numpy as np
import mdtraj as md
from numpy.fft import rfft, rfftfreq
from PMD import utils as utils
import os
import pandas as pd

def calc_com(traj_file,top_file):
    '''
       load traj without water and ions
       and return the com traj of each residue
    '''
    traj = md.load(traj_file,top=top_file)
    top = traj.topology
    N = traj.n_residues
    com_data = list()
    for i in range(N):
        atom_ids = top.select('backbone and resid '+str(i))
        temp_traj = traj.atom_slice(atom_ids)
        com_data.append(md.compute_center_of_mass(temp_traj))
        del temp_traj
    return com_data


def calc_fluc_spectra(com_data,time):
    '''
    Claculate the baseline of unbiased simulation
    Inputs:
        com_traj:   com trajectories used to calculate fluctuation power spectra
        time:       simulation time in ps
    Returns:
        xf:         FFT frequency used to plot figures
        fft_list:   list of FFT of com for each residue
    '''
    fft_list = list()
    for com_temp in com_data:
        com_temp1 = com_temp[:,0]
        com_temp2 = com_temp[:,1]
        com_temp3 = com_temp[:,2]
        yf1 = rfft(com_temp1-np.mean(com_temp1))
        yf2 = rfft(com_temp2-np.mean(com_temp2))
        yf3 = rfft(com_temp3-np.mean(com_temp3))
        yf = np.sqrt(yf1*yf1+yf2*yf2+yf3*yf3)
        fft_list.append(yf)
    N = len(com_data[0])
    xf = rfftfreq(N,time/N)
    return xf, fft_list
    


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
        atom_ids = utils.get_atom_ids(pdb_file,i+1,False)
        atom_ids_str_list = [str(id+1) for id in atom_ids]
        atom_mask = ','.join(atom_ids_str_list)
        out_file = 'res_'+str(i+1)+'_vac.dat'
        f.write('velocityautocorr @'+atom_mask+' out '+out_file+' tstep 0.01 norm\n')
    f.close()
    return None

def calc_dos(pdb_file,path):
    '''
    Input:
        path:   path to the folder where you store res_i_vac.dat
    Returns:
        xf:         dos frequency used to plot the ODS
        dos_list:   Density of states for each residue (use backbone atoms)
    '''
    res_num = utils.get_res_num(pdb_file)
    dos_list = list()
    for i in range(res_num):
        file_name_temp = path+'/res_'+str(i+1)+'_vac.dat'
        tab_temp = pd.read_table(file_name_temp,'\s+').values
        yf = rfft(tab[:,1]-np.mean(tab[:,1]))
        dos_list.append(yf)
        if i==0:
            N = len(tab_temp[:,0])
            time = tab_temp[1,0]+tab_temp[-1,0]
    xf = rfftfreq(N,time/N)
    return xf,dos_list

def pick_peak(xf,y0,y1,freq=1,ratio_threshold):
    '''
        pick the excited residue based on y1/y0 ration
        Returns the ratio list and excited residue list
    '''
    L = len(xf)
    idx = np.where(xf==freq)
    ratio = np.zeros(L)
    new_residues = list()
    for i in range(L):
        ratio[i] = y1[i][idx]/y0[i][idx]
        if ratio[i]>ratio_threshold:
            new_residues.append(i+1)
    return raio,new_residues
    