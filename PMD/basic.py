import numpy as np
import mdtraj as md
import os
from PMD import utils
from PMD import analysis as ana
import pickle
from prody import *
from PMD import nm
from itertools import cycle
import sys

class base:
    def __init__(self,pdb_file:str,path_dir:str,time:int,temperature:int,freq=1,force=100,init=False,N=10,velocity=False,GPU=[0,1,2,3]):
        '''
        construction function of base class
        pdb_file:       pdb_file name with same topology
        path_dir:       dir where the calculation is donw
        time:           simulation time
        temperature:    temperature of simulation (for powder simulation choose 120)
        freq:           pumped frequency. freq=1 means 1 ps per circle.   
        '''
        self.pdb = pdb_file
        self.n_res = utils.get_res_num(pdb_file)
        self.time = time
        self.temperature = temperature
        self.freq = freq
        self.force = force
        self.xf = None
        self.mean = None  
        self.std = None
        self.path_dir = path_dir
        if not os.path.exists(path_dir):
            os.system('mkdir '+path_dir)
        if init==False:
            self.index = 1
        else:
            self.index = 0
            self.control_simulation(N=N,velocity=velocity,GPU=GPU)
            self.index = 1
        
    
    def control_simulation(self,N=10,velocity=False,GPU=[0,1,2,3]):
        '''
        Run control simulation for N times
        And save the control_fps.pkl file for peak picking
        '''
        in_file = self.write_prod_input_file(None,velocity=velocity)
        f = open(self.path_dir+'/run_control.sh','w')
        f.write('#!/bin/bash\n')
        gpus = cycle(GPU)
        n_gpu = len(GPU)
        for i in range(N):
            out_file = self.path_dir+ '/control'+str(i)+'.out'
            x_file = self.path_dir + '/control'+str(i)+'.mdcrd'
            v_file = self.path_dir + '/control'+str(i)+'.mdvel'
            f.write('export CUDA_VISIBLE_DEVICES='+str(next(gpus))+'\n')
            if velocity==False:
                f.write('pmemd.cuda -O -i '+in_file+' -o '+out_file+' -p mol.prmtop -c equil.rst -x '+x_file)
            else:
                f.write('pmemd.cuda -O -i '+in_file+' -o '+out_file+' -p mol.prmtop -c equil.rst -x '+x_file+' -v '+v_file)
            if (i%n_gpu)!=(n_gpu-1) and i!=(N-1):
                f.write(' &')
            f.write('\n')
        f.close()
        os.system('bash '+self.path_dir+'/run_control.sh')
        file_name= self.path_dir+'/ctj_control.in'
        fps_list = list()
        for i in range(N):
            trajin = 'trajin '+self.path_dir+'/control'+str(i)+'.mdcrd'
            trajout = 'trajout '+self.path_dir+'/control'+str(i)+'_nw.xtc'
            f = open(file_name,'w')
            f.write('parm mol.prmtop\n')
            f.write(trajin+'\n')
            f.write('strip :WAT\n')
            f.write('strip :Na+\n')
            f.write('strip :Cl-\n')
            f.write(trajout+'\n')
            f.close()
            os.system('cpptraj -i '+file_name)
            traj_temp,rmsf_temp = ana.calc_com(self.path_dir+'/control'+str(i)+'_nw.xtc','mol_nw.prmtop',bb=False)
            xf,fps_temp = ana.calc_fps(traj_temp,rmsf_temp,time=self.time)
            fps_list.append(fps_temp)
        f = open(self.path_dir+'/control_fps.pkl','wb')
        pickle.dump([xf,fps_list],f)
        f.close()
        os.system('rm '+self.path_dir+'/*mdcrd')
        return None

    def write_plumed_file(self,atom_ids,vec,res_index=None):
        '''
        write plumed input file 
        return file name of plumed input file
        '''
        index = self.index
        if 3*len(atom_ids)!=len(vec):
            sys.exit('Atom numbers not match with vector length')
        if res_index is None:
            file_name = self.path_dir+'/pump'+str(index)+'.dat'
        else:
            file_name = self.path_dir+'/pump'+str(index)+'_'+res_index+'.dat'
        atoms = [str(atom_id+1) for atom_id in atom_ids]
        atoms_str = ','.join(atoms)
        k = str(round(self.force/1.661,2))
        omega = str(round(self.freq*6.28,2))
        x_weights_array = vec[0::3]
        x_weight = np.sum(x_weights_array)
        x_weight = str(x_weight.round(2))

        y_weights_array = vec[1::3]
        y_weight = np.sum(y_weights_array)
        y_weight = str(y_weight.round(2))

        z_weights_array = vec[2::3]
        z_weight = np.sum(z_weights_array)
        z_weight = str(z_weight.round(2))

        x_weights_list = [str(weight+0.001) for weight in x_weights_array]
        y_weights_list = [str(weight+0.001) for weight in y_weights_array]
        z_weights_list = [str(weight+0.001) for weight in z_weights_array]

        x_weights_str = ','.join(x_weights_list)
        y_weights_str = ','.join(y_weights_list)
        z_weights_str = ','.join(z_weights_list)

        f=open(file_name,'w')
        f.write('t: TIME\n')
        f.write('cx: CENTER ATOMS='+atoms_str+' WEIGHTS='+x_weights_str+' NOPBC\n')
        f.write('cy: CENTER ATOMS='+atoms_str+' WEIGHTS='+y_weights_str+' NOPBC\n')
        f.write('cz: CENTER ATOMS='+atoms_str+' WEIGHTS='+z_weights_str+' NOPBC\n')
        f.write('k: MATHEVAL ARG=t VAR=t FUNC='+k+'*'+'cos(t*'+omega+') PERIODIC=NO\n')
        f.write('Px: POSITION ATOM=cx\n')
        f.write('Py: POSITION ATOM=cy\n')
        f.write('Pz: POSITION ATOM=cz\n')
        f.write('V: MATHEVAL ARG=Px.x,Py.y,Pz.z,k VAR=x,y,z,k FUNC=('+x_weight+'*x+'+y_weight+'*y+'+z_weight+'*z)*k PERIODIC=NO\n')
        f.write('BV: BIASVALUE ARG=V\n')
        f.close()
        return file_name


    def write_prod_input_file(self,plumed_in,velocity=False):
        '''
        write prod input files
        return input file name
        '''
        index = self.index
        nsteps = str(int(1000*self.time))
        file_name = self.path_dir+'/prod'+str(index)+'.in'
        f = open(file_name,'w')
        f.write(file_name+'\n')
        f.write('&cntrl\n')
        f.write('  imin=0,irest=0,ntx=1,\n  nstlim='+nsteps+',dt=0.001,\n  ntc=2,ntf=2,\n')
        f.write('  cut=8.0, ntb=2, ntp=1, taup=2.0,\n')
        f.write('  ntpr=5000,ntwx=10,')
        if velocity:
            f.write('ntwv=10,')
        f.write('\n')
        f.write('  ntt=3,gamma_ln=2.0,temp0='+str(self.temperature)+',ig=-1')
        if index>0:
            f.write(',\n  plumed=1,plumedfile=\''+plumed_in+'\'\n/\n')
        else:
            f.write('\n/\n')
        f.close()
        return file_name


    def write_cpptraj_input_file(self):
        '''
        write cpptraj inpute file for extract waters and Na+ and CL-
        return the file name of cpptraj_file
        '''
        index=self.index
        file_name= self.path_dir+'/ctj'+str(index)
        trajin = 'trajin '+self.path_dir+'/prod'+str(index)+'.mdcrd'
        trajout = 'trajout '+self.path_dir+'/prod'+str(index)+'_nw.xtc'
        file_name=file_name+'.in'
        f = open(file_name,'w')
        f.write('parm mol.prmtop\n')
        f.write(trajin+'\n')
        f.write('strip :WAT\n')
        f.write('strip :Na+\n')
        f.write('strip :Cl-\n')
        f.write(trajout+'\n')
        f.close()
        return file_name

    def run(self,prod_in,cuda='0',velocity=False):
        index = self.index
        out_file = self.path_dir+ '/prod'+str(index)+'.out'
        x_file = self.path_dir + '/prod'+str(index)+'.mdcrd'
        v_file = self.path_dir + '/prod'+str(index)+'.mdvel'
        f = open(self.path_dir+'/run.sh','w')
        f.write('#!/bin/bash\n')
        f.write('export CUDA_VISIBLE_DEVICES='+cuda+'\n')
        if velocity==False:
            f.write('pmemd.cuda -O -i '+prod_in+' -o '+out_file+' -p mol.prmtop -c equil.rst -x '+x_file+'\n')
        else:
            f.write('pmemd.cuda -O -i '+prod_in+' -o '+out_file+' -p mol.prmtop -c equil.rst -x '+x_file+' -v '+v_file+'\n')
        f.close()
        os.system('bash '+self.path_dir+'/run.sh')
        return None

    def run_cpptraj(self):
        ctj_in = self.write_cpptraj_input_file()
        os.system('cpptraj -i '+ctj_in)
        return None

    def dump_traj(self,top_file='mol_nw.prmtop'):
        '''
        save traj information into pkl save time for further loading
        '''
        index = self.index
        traj_file = self.path_dir+'/prod'+str(index)+'_nw.xtc'
        f = open(self.path_dir+'/prod'+str(index)+'.pkl','wb')
        com_data,rmsf_data = ana.calc_com(traj_file,top_file)
        pickle.dump([com_data,rmsf_data],f)
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
