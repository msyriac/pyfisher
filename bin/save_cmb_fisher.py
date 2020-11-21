from __future__ import print_function
from orphics import maps,io,cosmology,stats,mpi
import numpy as np
import os,sys
import pyfisher

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Do a thing.')
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-i","--input",type=str,help="Input root",required=True)
required_args.add_argument("-p","--param-file",type=str,help="Parameter file",required=True)
args = parser.parse_args()

out_name = f"{args.input}"
param_dat = np.genfromtxt(args.param_file,dtype=None,encoding='utf-8',delimiter=',')
_,fids = pyfisher.get_jobs(args.param_file,exclude=None)


"""
Useful Fisher matrices
Planck TT fsky=0.65 for 2 < ell < 30
Planck TT,TE,EE fsky=0.6 - fsky_exp for 30 < ell < 2500
A tau prior


"""


# Planck noise
lmax = 3000
ells = np.arange(lmax)
beams_T =  [33.,23.,14.,10.,7.,5.,5.]
uK_arcmins_T = [145.,149.,137.,65.,43.,66.,200.]
beams_P =  [14.,10.,7.,5.,5.]
uK_arcmins_P = [450.,103.,81.,134.,406.]
Ns_TT = np.asarray([(uK_arcmin*np.pi/180./60.)**2./maps.gauss_beam(ells,fwhm)**2. for uK_arcmin,fwhm in zip(uK_arcmins_T,beams_T)])
Ns_PP = np.asarray([(uK_arcmin*np.pi/180./60.)**2./maps.gauss_beam(ells,fwhm)**2. for uK_arcmin,fwhm in zip(uK_arcmins_P,beams_P)])
N_TT = 1./(1./Ns_TT).sum(axis=0)
N_PP = 1./(1./Ns_PP).sum(axis=0)
nls = {}
nls['TT'] = maps.interp(ells,N_TT)
nls['EE'] = maps.interp(ells,N_PP)
nls['BB'] = maps.interp(ells,N_PP)

# CMB Cls

def load_theory_dict(fname):
    cls = {}
    ells,tt,ee,bb,te,kk = np.loadtxt(fname,unpack=True)
    cls['TT'] = maps.interp(ells,tt)
    cls['EE'] = maps.interp(ells,ee)
    cls['BB'] = maps.interp(ells,bb)
    cls['TE'] = maps.interp(ells,te)
    cls['kk'] = maps.interp(ells,kk)
    return cls

cls = load_theory_dict(f'{args.input}_fiducial.txt')

param_list = list(fids.keys())
dcls = {}
for param in param_list:
    dcls[param] = load_theory_dict(f'{args.input}_{param}_deriv.txt')



fsky1 = 0.65
fsky2 = 0.65

bin_edges = np.arange(2,31)
specs = ['TT']
F1 = pyfisher.band_fisher(param_list,bin_edges,specs,cls,nls,dcls) 
bin_edges = np.arange(30,2501)
specs = ['TT','EE','TE']
F2 = pyfisher.band_fisher(param_list,bin_edges,specs,cls,nls,dcls)

F = fsky1*F1 + fsky2*F2
F.delete(['nnu','mnu','ok','w0','wa'])
F.add_prior('tau',0.01)
print(F)
print(F.sigmas())

