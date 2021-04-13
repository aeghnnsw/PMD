import pickle
import numpy as np
import mdtraj as md
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

def get_atom_ids(pdb_file,res_ids,calpha=True):
    '''
    Takes pdb file name, and res_ids as input and returns the atoms ids
    Request: prepare the pdb that has added H atoms
    Inputs:
        pdb_file:   string
        res_ids:    list of integers
        calpha:     bool, default to be True, when true, only select calpha atoms
    Return:
        atom_ids:   list of integers
    '''
    pdb = md.load_pdb(pdb_file)
    top = pdb.topology
    atom_ids = list()
    if isinstance(res_ids,int):
        if calpha:
            atom_ids=top.select('name CA and resid '+str(res_ids-1))
        else:
            atom_ids=top.select('backbone and resid '+str(res_ids-1))
    else:
        for res_id in res_ids:
            if calpha:
                atom_ids.extend(top.select('name CA and resid '+str(res_id-1)))
            else:
                atom_ids.extend(top.select('backbone and resid '+str(res_id-1)))
    return atom_ids

def get_res_num(pdb_file):
    pdb = md.load_pdb(pdb_file)
    return pdb.n_residues

def get_dist(pdb_file,res1,res2):
    pdb = md.load_pdb(pdb_file)
    top = pdb.topology
    id1 = top.select('name CA and resid '+str(res1-1))
    id2 = top.select('name CA and resid '+str(res2-1))
    x1 = pdb.xyz[0][id1]
    x2 = pdb.xyz[0][id2]
    delta = x1-x2
    return np.sqrt(np.sum(delta**2))
    
def get_AM(cm='connection_map.pkl',pdb='mol.pdb'):
    f = open(cm,'rb')
    m = pickle.load(f)
    f.close()
    L = get_res_num(pdb)
    AM = np.zeros([L,L])
    for i in m:
        for j in m[i]:
            AM[i-1][j-1]=1
    for i in range(L):
        AM[i][i] = 1
    for i in range(1,L):
        AM[i][i-1]=1
        AM[i-1][i]=1
    return AM

def calc_shortest_path(AM):
    graph = csr_matrix(AM)
    dist,pred = shortest_path(graph,return_predecessors=True)
    return dist,pred

def get_shortest_path(res_i,res_j,pred):
    path_temp = [res_j]
    i = res_i-1
    j = res_j-1
    temp = j
    while temp!=i:
        temp = pred[i,temp]
        path_temp.append(temp+1)
    path = list()
    while len(path_temp)!=0:
        path.append(path_temp.pop())
    return path

    
    









    
    
