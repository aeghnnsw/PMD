import numpy as np
import mdtraj as md
from numpy.fft import rfft, rfftfreq
from PMD import utils as utils
import os
import pandas as pd

def calc_com(traj_file,top_file,bb=True):
    '''
       load traj without water and ions
       and return the com traj of each residue
    '''
    traj = md.load(traj_file,top=top_file)
    top = traj.topology
    N = traj.n_residues
    com_data = list()
    for i in range(N):
        if bb==True:
            atom_ids = top.select('backbone and resid '+str(i))
        else:
            atom_ids = top.select('name CA and resid '+str(i))
        temp_traj = traj.atom_slice(atom_ids)
        com_data.append(md.compute_center_of_mass(temp_traj))
        del temp_traj
    return com_data


def calc_fps(traj_data,time):
    '''
    Claculate the baseline of unbiased simulation
    Inputs:
        traj:       trajectories used to calculate fluctuation power spectra, list
        time:       simulation time in ps
    Returns:
        xf:         FFT frequency used to plot figures
        fps_list:   fluctuation power spectra list for each residue
    '''
    fps_list = list()
    for traj_temp in traj_data:
        traj_temp1 = traj_temp[:,0]
        traj_temp2 = traj_temp[:,1]
        traj_temp3 = traj_temp[:,2]
        yf1 = np.abs(rfft(traj_temp1-np.mean(traj_temp1)))
        yf2 = np.abs(rfft(traj_temp2-np.mean(traj_temp2)))
        yf3 = np.abs(rfft(traj_temp3-np.mean(traj_temp3)))
        yf = yf1*yf1+yf2*yf2+yf3*yf3
        fps_list.append(yf)
    N = len(traj_data[0])
    xf = rfftfreq(N,time/N)
    return xf, fps_list
    


def write_cpptraj_vac_input_file(pdb_file,crd_file,vel_file,top_file,path):
    '''
    write cpptraj input file to calculate VAC for each residue
    Input:
        crd_file:   crd file name
        vel_file;   vel file name
        top_file:   topology file name
        path:       path directory
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
        out_file = path+'/res_'+str(i+1)+'_vac.dat'
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
        yf = np.abs(rfft(tab_temp[:,1]-np.mean(tab_temp[:,1])))
        dos_list.append(yf)
        if i==0:
            N = len(tab_temp[:,0])
            time = tab_temp[1,0]+tab_temp[-1,0]
    xf = rfftfreq(N,time/N)
    return xf,dos_list

def pick_peak(xf,y0,y1,ratio_threshold=2,top=5,window=1,freq=1):
    '''
        pick the excited residue based on y1/y0 ratio
        And the difference between y1-y0, (add somw window to alleviate problems raised by small y0 values)
        Returns the excited residue list
    '''
    L = len(y0)
    idx = np.where(xf==freq)[0][0]
    dif = np.zeros(L)
    new_residues = list()
    for i in range(L):
        dif_temp = y1[i][idx]-y0[i][idx]
        if window>0:
            for j in range(window):
                dif_temp = dif_temp + (y1[i][idx+j+1] - y0[i][idx+j+1])
                dif_temp = dif_temp + (y1[i][idx-j-1] - y0[i][idx-j-1])
        dif[i] = dif_temp
    sort_dif = np.sort(dif)
    dif_threshold = sort_dif[-top]
    new_idx = np.where(dif>=dif_threshold)
    for j in new_idx[0]:
        if (y1[j][idx]/y0[j][idx])>ratio_threshold:
            new_residues.append(j+1)
    return new_residues
    
