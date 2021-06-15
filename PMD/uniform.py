import numpy as np
import mdtraj as md
import os
from prody.dynamics import mode
from PMD import utils
from PMD import analysis as ana
from PMD.basic import base
from prody import *
from PMD import nm
from itertools import cycle
import sys

class unif(base):
    def __init__(self,pdb_file:str,path_dir:str,time:int,temperature:int,vec:list,freq=1,force=100,init=False,N=10,velocity=False,GPU=[0,1,2,3]):
        super(unif,self).__init__(pdb_file=pdb_file,path_dir=path_dir,time=time,temperature=temperature,freq=freq,force=force,init=init,N=N,velocity=velocity,GPU=GPU)
        if len(vec)!=3:
            sys.exit('please input a vector with three dimension')
        vec = np.array(vec)
        vec = vec/np.sqrt(np.sum(vec*vec))
        self.vec = vec
        

    def calc_pump_vec(self,atom_ids):
        N = len(atom_ids)
        pump_vec = np.zeros(N*3)
        pump_vec[0::3] = self.vec[0]
        pump_vec[1::3] = self.vec[1]
        pump_vec[2::3] = self.vec[2]
        return pump_vec.round(2)

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
        vec = self.calc_pump_vec(atom_ids)
        plumed_in = self.write_plumed_file(atom_ids,vec,res_index=res_index)
        prod_in = self.write_prod_input_file(plumed_in)
        self.run(prod_in,cuda=cuda)
        self.run_cpptraj()
        self.dump_traj()
        if rm_file:
            os.system('rm '+self.path_dir+'/*mdcrd')
        self.index+=1
        return None
