#   The main program to run loop

import PMD.utils as utils
import os
import argparse


parser = argparse.ArgumentParser(description="The main process to run PMD\n Loop to find the pathway between residue i to residue j")
parser.add_argument('--pdb',action='store',type=str,required=True,help='The pdb file path to be used')
parser.add_argument('--start_res',action='store',type=str,required=True,help='The ids of starting residues, saperate by \',\' no space, for example 1,2,3')
parser.add_argument('--end_res',action='store',type=str,required=True,help='The ids of end residues, saperate by \',\' no space')
parser.add_argument('--force',action='store',type=float,required=True,help='The amplitude of oscilation force, unit pN')
parser.add_argument('--frequency',action='store',type=float,required=True,help='Frequency of oscilation force, unit 1ps-1 or 10^12 Hz')

args = parser.parse_args()
start_res = args.start_res.split(',')
end_res = args.end_res.split(',')
start_res = [int(id) for id in start_res]
end_res = [int(id) for id in end_res]

index = 0
start_ids = utils.get_atom_ids(start_res)
sys.write_plumed_files(start_ids,args.frequency,args.force,index)







