from __future__ import print_function
from orphics import maps,io,cosmology
from enlib import enmap
import numpy as np
import os,sys


broot = lambda param : "output/BAO_highAcc_DESI2_szar_step_0.01_dfk_%s.csv" % param
broot2 = lambda param : "output/BAO_highAcc_DESI2_szar__dfk_%s.csv" % param

params = "H0,ombh2,omch2,tau,As,ns,mnu,w,wa".split(',')
zs = np.array([.15,.25,.35,.45,.55,.65,.75,.85,.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85])
for param in params:

    pl = io.Plotter()
    d = np.loadtxt(broot(param))
    d2 = np.loadtxt(broot2(param))
    pl.add(zs,d)
    pl.add(zs,d2,ls="--")
    pl.done(io.dout_dir+"dbao_%s.png" % param)
