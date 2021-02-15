import numpy as np
import mdtraj as md

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

    