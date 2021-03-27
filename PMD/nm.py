import numpy as np
import mdtraj as md
import os
from prody import *
from PMD import utils as utils
from PMD import analysis as ana
import pickle

class nm:
    def __init__(self,pdb_file,path_dir):
        self.pdb = pdb_file
        self.n_res = utils.get_res_num(pdb_file)
        self.index = 0
        self.xf = None
        self.y0 = None        
        self.path_dir = path_dir
        if not os.path.exists(path_dir):
            os.system('mkdir '+path_dir)

    def calc_nm_vec(self,res_ids,mode_number=0):
        '''
        Inputs:
            pdb_file:       string
            res_ids:        list of integers or integer
            mode_number:    Mode number indicate which normal mode to use, defautl=0
        Returns:
            nm_vec:     Normal modes vector for input res_ids, shape: n_res*3
        '''
        pdb_file = self.pdb
        if isinstance(res_ids,int):
            id_str = str(res_ids)
        else:
            id_str_list = [str(i) for i in res_ids]
            id_str = ' '.join(id_str_list)
        protein = parsePDB(pdb_file)
        calphas = protein.select('calpha')
        target_id = protein.select('calpha and resnum '+id_str).getResnums()
        N = len(target_id)
        atom_list = list(calphas.getResnums())
        tar_res_id = [atom_list.index(i) for i in target_id]
        vec_ids = list()
        for temp_id in tar_res_id:
            vec_ids.append(temp_id*3)
            vec_ids.append(temp_id*3+1)
            vec_ids.append(temp_id*3+2)
        anm = ANM('ANM analysis')
        anm.buildHessian(calphas)
        anm.calcModes()
        mode = anm[mode_number]
        eigvec = mode.getEigvec()
        vec = eigvec[vec_ids]
        vec = vec/np.sqrt(np.sum(vec*vec))
        vec = vec*np.sqrt(N)
        return vec.round(2)


    def calc_nm_fluc(self,mode_number=0):
        '''
        calculate the fluctuations for each residue at given mode number
        Returns:
            The fluctuation vec (length is the number of residues)
        '''
        pdb_file = self.pdb
        protein = parsePDB(pdb_file)
        calphas = protein.select('calpha')
        anm = ANM('ANM analysis')
        anm.buildHessian(calphas)
        anm.calcModes()
        mode = anm[mode_number]
        eigvec = mode.getEigvec()
        disp = eigvec.reshape(-1,3)
        fluc = np.sqrt(np.sum(disp*disp,1))
        return fluc


    def write_plumed_file(self,atom_ids,vec,frequency=1,force=100,res_index=None):
        '''
        write plumed input files for normal modes input 
        Inputes:
            vec:    the eigenvec of the normal modes for the selected atoms
        '''
        index = self.index
        if res_index is None:
            file_name = self.path_dir+'/pump'+str(index)+'_nm.dat'
        else:
            file_name = self.path_dir+'/pump'+str(index)+'_'+res_index+'nm.dat'
        atoms = [str(atom_id+1) for atom_id in atom_ids]
        atoms_str = ','.join(atoms)
        k = str(round(force/1.661,2))
        omega = str(round(frequency*6.28,2))

        x_weights_array = vec[0::3]
        x_weight = np.sqrt(np.sum(x_weights_array*x_weights_array))
        x_weight = str(x_weight.round(2))
        y_weights_array = vec[1::3]
        y_weight = np.sqrt(np.sum(y_weights_array*y_weights_array))
        y_weight = str(y_weight.round(2))
        z_weights_array = vec[2::3]
        z_weight = np.sqrt(np.sum(z_weights_array*z_weights_array))
        z_weight = str(z_weight.round(2))
    
        x_weights_list = [str(weight) for weight in x_weights_array]
        y_weights_list = [str(weight) for weight in y_weights_array]
        z_weights_list = [str(weight) for weight in z_weights_array]

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


    def write_prod_input_file(self,time,plumed_in,velocity=False):
        '''
        write prod input files
        Inputs:
            time:   int in ps
            index:  if index=0 it is unbiased simulation
        return input file name
        '''
        index = self.index
        nsteps = str(int(1000*time))
        file_name = self.path_dir+'/prod'+str(index)+'_nm.in'
        f = open(file_name,'w')
        f.write(file_name+'\n')
        f.write('&cntrl\n')
        f.write('  imin=0,irest=1,ntx=5,\n  nstlim='+nsteps+',dt=0.001,\n  ntc=2,ntf=2,\n')
        f.write('  cut=8.0, ntb=2, ntp=1, taup=2.0,\n')
        f.write('  ntpr=5000,ntwx=10,')
        if velocity:
            f.write('ntwv=10,')
        f.write('\n')
        f.write('  ntt=3,gamma_ln=2.0,temp0=300.0,ig=-1')
        if index>0:
            f.write(',\n  plumed=1,plumedfile=\''+plumed_in+'\'\n/\n')
        else:
            f.write('\n/\n')
        f.close()
        return file_name

    def write_cpptraj_input_file(self):
        '''
        write cpptraj inpute file for extract waters and Na+ and CL-
        return the file name
        '''
        index=self.index
        file_name= self.path_dir+'/ctj'+str(index)
        trajin = 'trajin '+self.path_dir+'/prod'+str(index)+'_nm.mdcrd'
        trajout = 'trajout '+self.path_dir+'/prod'+str(index)+'_nm_nw.xtc'
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

    def run_cpptraj(self):
        '''
        run cpptraj to remove waters and ions of a trajectory and save the file
        Inputs:
            ctj_in: str name of cpptraj input file
        '''
        ctj_in = self.write_cpptraj_input_file()
        os.system('cpptraj -i '+ctj_in)
        return None

    def strip_topology_wat(self):
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

    def run(self,in_file,cuda='0',velocity=False):
        index = self.index
        out_file = self.path_dir+ '/prod'+str(index)+'_nm.out'
        x_file = self.path_dir + '/prod'+str(index)+'_nm.mdcrd'
        v_file = self.path_dir + '/prod'+str(index)+'_nm.mdvel'
        f = open(self.path_dir+'/run.sh','w')
        f.write('#!/bin/bash\n')
        f.write('export CUDA_VISIBLE_DEVICES='+cuda+'\n')
        if velocity==False:
            f.write('pmemd.cuda -O -i '+in_file+' -o '+out_file+' -p mol.prmtop -c equil.rst -x '+x_file+'\n')
        else:
            f.write('pmemd.cuda -O -i '+in_file+' -o '+out_file+' -p mol.prmtop -c equil.rst -x '+x_file+' -v '+v_file+'\n')
        f.close()
        os.system('bash '+self.path_dir+'/run.sh')

        
    def pump(self,res_ids,mode=0,frequency=1,force=100,time=500,velocity=False,cuda='0',ratio_threshold=2,top=5,window=1,rm_file=True):
        '''
        Run pumped md and return the pumped residues if index>0
        '''
        pdb_file = self.pdb
        index = self.index
        if index==0:
            if os.path.exists('control.pkl'):
                f = open('control.pkl','rb')
                self.xf,self.y0 = pickle.load(f)
                f.close()
                self.index += 1
            else:
                in_file = self.write_prod_input_file(time,None,velocity)
                self.run(in_file,cuda,velocity=velocity)
                self.run_cpptraj()
                self.strip_topology_wat()
                traj_file = self.path_dir+'/prod'+str(index)+'_nm_nw.xtc'
                top_file = 'mol_nw.prmtop'
                traj_data = ana.calc_com(traj_file,top_file)
                self.xf,self.y0 = ana.calc_fps(traj_data,time)
                f = open('control.pkl','wb')
                pickle.dump([self.xf,self.y0],f)
                f.close()   
                self.index +=1
        index = self.index
        atom_ids = utils.get_atom_ids(pdb_file,res_ids)
        vec = self.calc_nm_vec(res_ids,mode)
        plumed_in = self.write_plumed_file(atom_ids,vec,frequency,force)
        in_file = self.write_prod_input_file(time,plumed_in,velocity)
        self.run(in_file,cuda,velocity=velocity)
        self.run_cpptraj()
        traj_file = self.path_dir+'/prod'+str(index)+'_nm_nw.xtc'
        top_file = 'mol_nw.prmtop'
        traj_data = ana.calc_com(traj_file,top_file)
        _,y_temp = ana.calc_fps(traj_data,time)
        nr = ana.pick_peak(self.xf,self.y0,y_temp,ratio_threshold=ratio_threshold,top=top,window=window,freq=frequency)
        self.index +=1
        if rm_file:
            os.system('rm '+self.path_dir+'/*mdcrd')
        return nr

    def pump_at_res(self,res_ids,mode=0,frequency=1,force=100,time=500,velocity=False,cuda='0',ratio_threshold=2,top=5,window=1,rm_file=True):
        '''
        Run pumped md and return the pumped residues if index>0
        '''
        pdb_file = self.pdb
        self.index = -1
        if isinstance(res_id,int):
            res_index = str(res_ids)
        if isinstance(res_id,list):
            res_index = str()
            for temp_id in res_ids:
                res_index = res_index+str(temp_id)+'_'
        atom_ids = utils.get_atom_ids(pdb_file,res_ids)
        vec = self.calc_nm_vec(res_id,mode)
        plumed_in = self.write_plumed_file(atom_ids,vec,frequency=frequency,force=force,res_index=res_index)
        in_file = self.write_prod_input_file(time,plumed_in,velocity)
        self.run(in_file,cuda,velocity=velocity)
        self.run_cpptraj()
        traj_file = self.path_dir+'/prod'+str(index)+'_nm_nw.xtc'
        top_file = 'mol_nw.prmtop'
        traj_data = ana.calc_com(traj_file,top_file)
        _,y_temp = ana.calc_fps(traj_data,time)
        nr = ana.pick_peak(self.xf,self.y0,y_temp,ratio_threshold=ratio_threshold,top=top,window=window,freq=frequency)
        if rm_file:
            os.system('rm '+self.path_dir+'/*mdcrd')
        return nr
        

        

