from __future__ import print_function
import numpy as np
import os,sys
import pyfisher
from pyfisher import mpi
import argparse

def parse_args(args):
    # Parse command line
    parser = argparse.ArgumentParser(description='Do a thing.')
    parser.add_argument("--exp", type=str,default='planck',help='Experiment name.')
    parser.add_argument("--Lmin",     type=int,  default=4,help="A description.")
    parser.add_argument("--Lmax",     type=int,  default=400,help="A description.")
    parser.add_argument("--fsky",     type=float,  default=0.65,help="A description.")
    args = parser.parse_args()
    return args

def main():
    args = parse_args(sys.argv[1:])
    run_lensing_demo(args.Lmin,args.Lmax,args.exp,args.fsky,test=False)

def run_lensing_demo(Lmin,Lmax,exp,fsky,test=False):

    # Load default fiducial parameters
    fids = pyfisher.get_fiducials()

    # Load a pre-calculated BOSS BAO Fisher
    bao = pyfisher.get_saved_fisher('boss_bao')

    # Load a pre-calculated CMB lensing noise curve
    ells,nls = pyfisher.get_lensing_nl(exp)
    # Calculate a CMB lensing Fisher
    bin_edges = np.arange(Lmin,Lmax)
    lens = pyfisher.get_lensing_fisher(bin_edges,ells,nls,fsky)

    # Planck lens + BAO (s8, om, H0 parameterization)
    F = lens+bao
    F.delete(['w0','wa','ok','nnu','tau','mnu'])
    s8 = pyfisher.get_s8(zs=[0.],params=fids)[0]
    fids['s8'] = s8 
    fids['om'] = (fids['omch2'] + fids['ombh2'])/(fids['H0']/100)**2.
    F1 = pyfisher.reparameterize(F,['om','s8','H0','ns','ombh2'],fids,verbose=False)
    F1.add_prior('ns',0.02)
    F1.add_prior('ombh2',0.0005)
    sigmas = F1.sigmas()
    print("Planck lens + BAO (s8, om, H0 parameterization)")
    for p in ['s8','om','H0']:
        print(f'{p} = {fids[p]:.03f}+-{sigmas[p]:.03f}')
    if test: 
        assert np.isclose(sigmas['s8'],0.01539,rtol=1e-3)
        assert np.isclose(sigmas['om'],0.01903,rtol=1e-3)

    # Planck lens alone (s8om^0.25, H0 parameterization)
    fids['s8om0.25'] = fids['s8'] * fids['om']**0.25
    F = pyfisher.reparameterize(lens,['s8om0.25','H0','ns','ombh2'],fids,verbose=False)
    F.add_prior('ns',0.02)
    F.add_prior('ombh2',0.0005)
    sigmas = F.sigmas()
    print("Planck lens  (s8om^0.25, H0 parameterization)")
    for p in ['s8om0.25']:
        print(f'{p} = {fids[p]:.03f}+-{sigmas[p]:.03f}')
    if test:
        assert np.isclose(sigmas['s8om0.25'],0.01516,rtol=1e-3)


    # Planck lens + BAO + CMB (mnu)
    lcmb = pyfisher.get_saved_fisher('planck_lowell',0.65)
    hcmb = pyfisher.get_saved_fisher('planck_highell',0.65)
    F = lens + bao + lcmb + hcmb
    F.delete(['w0','wa','ok','nnu'])
    F.add_prior('tau',0.007)
    sigmas = F.sigmas()
    print("Planck lens + BAO + CMB (mnu)")
    for p in ['mnu']:
        print(f'{p} = {fids[p]:.03f}+-{sigmas[p]:.03f}')
    if test:
        assert np.isclose(sigmas['mnu'],0.0743,rtol=1e-3)
    else:
        pyfisher.contour_plot(F,fids,'contour.png',name=None)

def test_lensing_demo():
    run_lensing_demo(4,400,'planck',0.65,test=True)

if __name__ == "__main__":
    main()

