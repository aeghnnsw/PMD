import numpy as np
import mdtraj as md
import os

from prody.dynamics import mode
from PMD import utils
from PMD import analysis as ana
from PMD.basic import base
from prody import *
from PMD import nm


class nm(base):
    def __init__(self,pdb_file:str,path_dir:str,time:int,temperature:int,mode_index=0,freq=1,force=100,init=False,N=10,velocity=False,GPU=[0,1,2,3]):
        super(nm,self).__init__(pdb_file=pdb_file,path_dir=path_dir,time=time,temperature=temperature,freq=freq,force=force,init=init,N=N,velocity=velocity,GPU=GPU)
        self.mode = mode_index
    
    def calc_nm_vec(self,res_ids):
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
        mode = anm[self.mode]
        eigvec = mode.getEigvec()
        vec = eigvec[vec_ids]
        vec = vec/np.sqrt(np.sum(vec*vec))
        vec = vec*np.sqrt(N)
        return vec.round(2)

    def pump_at_res(self,res_ids,cuda='0',rm_file=True):
        '''
        Run pumped MD given res_id
        '''
        pdb_file = self.pdb
        if isinstance(res_ids,int):
            res_index = str(res_ids)
        if isinstance(res_ids,list):
            res_index = str()
            for temp_id in res_ids:
                res_index = res_index+'_'+str(temp_id)
        atom_ids = utils.get_atom_ids(pdb_file,res_ids)
        vec = self.calc_nm_vec(res_ids)
        plumed_in = self.write_plumed_file(atom_ids,vec,res_index=res_index)
        prod_in = self.write_prod_input_file(plumed_in)
        self.run(prod_in,cuda=cuda)
        self.run_cpptraj()
        self.dump_traj()
        if rm_file:
            os.system('rm '+self.path_dir+'/*mdcrd')
        self.index+=1
        return None


