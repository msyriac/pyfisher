from __future__ import print_function
import numpy as np
import os,sys,shutil
import pyfisher
from pyfisher import mpi

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Do a thing.')
parser.add_argument("--lmax",     type=int,  default=3000,help="A description.")
parser.add_argument("--accurate", action='store_true',help='A flag.')
parser.add_argument("--linear", action='store_true',help='A flag.')
parser.add_argument("--lensed", action='store_true',help='A flag.')
parser.add_argument("--exclude-in-lensing",     type=str,  default='tau,r',help="Relative path to directory with experiment info.")
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-o","--output",type=str,help="Output root",required=True)
required_args.add_argument("-p","--param-file",type=str,help="Parameter file",required=True)
args = parser.parse_args()

out_name = pyfisher.prepare_output(args,"save_cmb.py CMB derivatives run")

lens_exclude = args.exclude_in_lensing.split(',')
shutil.copyfile(args.param_file,args.output+"/params.txt")
jobs,fids = pyfisher.get_param_info(args.param_file,exclude=None)
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

    theory = pyfisher.get_cls(params=pparams,lmax=args.lmax,accurate=args.accurate,engine='camb',de='ppf',nonlinear=not(args.linear)) 
    ells = np.arange(args.lmax)
    cfunc = theory.lCl if args.lensed else theory.uCl
    saves = []
    saves.append(ells)
    for spec in ['TT','EE','BB','TE']:
        saves.append( cfunc(spec,ells) )
    saves.append( theory.gCl('kk',ells) )

    if param is None:
        fname = f'{out_name}_cmb_fiducial.txt'
    else:
        fname = f'{out_name}_cmb_{param}_{ptype}.txt'
    hstr = ','.join([param,str(val)]) if param is not None else ""
    io.save_cols(fname,saves,header=hstr)

comm.Barrier()
if rank==0:
    def read(param,ud):
        filename = f'{out_name}_cmb_{param}_{ud}.txt'
        with open(filename,'r') as f:
            header = f.readline().strip()
            assert header[0]=="#"
            oparam,val = header[1:].split(',')
        data = np.loadtxt(filename)
        return oparam.strip(),float(val),data

    fiducial_theory = np.loadtxt(f'{out_name}_cmb_fiducial.txt')
    deriv_theory = {}
    for param in fids.keys():
        uparam,uval,udata = read(param,'u')
        dparam,dval,ddata = read(param,'d')
        assert param==uparam==dparam
        deriv = (udata[:,1:]-ddata[:,1:]) / (uval-dval)
        odat = udata.copy()
        odat[:,1:] = deriv.copy()
        if param in lens_exclude:
            print(f"Zeroing lensing derivative for {param}")
            odat[:,-1]  = 0
        fname = f'{out_name}_cmb_{param}_deriv.txt'
        np.savetxt(fname,odat)
