from __future__ import print_function
import numpy as np
import os,sys,shutil
import pyfisher
from pyfisher import mpi

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Do a thing.')
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-i","--input",type=str,help="Input root",required=True)
required_args.add_argument("-o","--output",type=str,help="Output root",required=True)
args = parser.parse_args()

out_name = pyfisher.prepare_output(args,"save_cmb_fisher.py Planck CMB Fishers run")

param_file = f'{args.input}/params.txt'
_,fids = pyfisher.get_param_info(param_file,exclude=None)

param_list = list(fids.keys())
F1 = pyfisher.get_planck_cmb_fisher(param_list,np.arange(2,31),['TT'],args.input,fsky=1)
F2 = pyfisher.get_planck_cmb_fisher(param_list,np.arange(30,2501),['TT','EE','TE'],args.input,fsky=1)

pyfisher.write_fisher(f'{out_name}_planck_low_ell_TT_fullsky.txt',F1,delim=',')
pyfisher.write_fisher(f'{out_name}_planck_high_ell_TTEETE_fullsky.txt',F2,delim=',')
