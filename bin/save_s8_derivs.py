from __future__ import print_function
from orphics import maps,io,cosmology,stats,mpi
import numpy as np
import os,sys,shutil
import pyfisher

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Do a thing.')
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-n","--nsteps",type=int,help="Number of steps",required=True)
required_args.add_argument("-o","--output",type=str,help="Output root",required=True)
required_args.add_argument("-p","--param-file",type=str,help="Parameter file",required=True)
parser.add_argument("--exclude",     type=str,  default='tau,r',help="Relative path to directory with experiment info.")
args = parser.parse_args()

output_root = pyfisher.prepare_output(args,"save_s8_derivs.py s8 deriv run")
exclude = args.exclude.split(',')
jobs,fids = pyfisher.get_param_info(args.param_file,exclude,get_range=True)
shutil.copyfile(args.param_file,args.output+"/params.txt")
njobs = len(jobs)
comm,rank,my_tasks = mpi.distribute(njobs)

for task in my_tasks:
    param,pmin,pmax = jobs[task]
    print(param,pmin,pmax)

    vals = np.linspace(pmin,pmax,args.nsteps)
    s8s = []
    pl = io.Plotter(xyscale='linlin',xlabel=pyfisher.latex_mapping[param],ylabel=pyfisher.latex_mapping['s8'])
    for val in vals:
        pparams = dict(fids)
        pparams[param] = val
        s8 = pyfisher.get_s8(zs=[0.],params=pparams)[0]
        s8s.append(s8)
    pl.add(vals,s8s,marker='o')
    p1,p0 = np.polyfit(vals, s8s, 1)
    pl.add(vals,p0+p1*vals,ls='--')
    pl.done(f'{output_root}s8s_{param}.png')
    np.savetxt(f'{output_root}_fitderiv_s8_wrt_{param}.txt')
