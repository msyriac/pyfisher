from __future__ import print_function
import numpy as np
import os,sys,shutil
import pyfisher
from pyfisher import mpi
import hmvec

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Do a thing.')
parser.add_argument("exp_name", type=str,help='Positional arg.')
parser.add_argument("--boss-include",     type=str,  default='6df,mgs,lowz,cmass',help="A description.")
parser.add_argument("--input-path",     type=str,  default='input',help="Relative path to directory with experiment info.")
parser.add_argument("--exclude",     type=str,  default='As,ns,tau,r,S8,sigma8',help="Relative path to directory with experiment info.")
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-o","--output",type=str,help="Output root",required=True)
required_args.add_argument("-p","--param-file",type=str,help="Parameter file",required=True)
args = parser.parse_args()


output_root = pyfisher.prepare_output(args,"save_bao.py BAO Fisher Matrix run")

exclude = args.exclude.split(',')
zs,sig_pers = pyfisher.load_bao_experiment_rs_dV_diagonal(args.exp_name,args.input_path,boss_include=args.boss_include.split(','))
jobs,fids = pyfisher.get_param_info(args.param_file,exclude)
shutil.copyfile(args.param_file,args.output+"/params.txt")
njobs = len(jobs)
comm,rank,my_tasks = mpi.distribute(njobs)

for task in my_tasks:
    param,val,ptype = jobs[task]
    print(param,val,ptype)
    pparams = dict(fids)
    if param is None:
        assert val is None
        assert ptype=='f'
    else:
        pparams[param] = val

    hmvec.Cosmology(pparams,engine='class')
    retval = hm.get_bao_rs_dV(zs)

    if param is None:
        fname = f'{output_root}_bao_fiducial.txt'
    else:
        fname = f'{output_root}_bao_{param}_{ptype}.txt'
    hstr = ','.join([param,str(val)]) if param is not None else ""
    np.savetxt(fname,retval,header=hstr)

comm.Barrier()
if rank==0:
    def read(param,ud):
        filename = f'{output_root}_bao_{param}_{ud}.txt'
        with open(filename,'r') as f:
            header = f.readline().strip()
            assert header[0]=="#"
            oparam,val = header[1:].split(',')
        data = np.loadtxt(filename)
        return oparam.strip(),float(val),data

    fiducial_theory = np.loadtxt(f'{output_root}_bao_fiducial.txt')
    deriv_theory = {}
    for param in fids.keys():
        uparam,uval,udata = read(param,'u')
        dparam,dval,ddata = read(param,'d')
        assert param==uparam==dparam
        deriv_theory[param] = (udata-ddata) / (uval-dval)
    Fmat = pyfisher.get_bao_fisher_rs_dV_diagonal(list(fids.keys()),deriv_theory,fiducial_theory,sig_pers)
    pyfisher.write_fisher(f'{output_root}_bao_fisher.txt',Fmat,delim=',')
