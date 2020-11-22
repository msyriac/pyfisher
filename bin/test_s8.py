from __future__ import print_function
from orphics import maps,io,cosmology,stats,lensing
from pixell import enmap
import numpy as np
import os,sys
import pyfisher

fids = pyfisher.get_fiducials()
bao = pyfisher.get_saved_fisher('desi_bao')

bin_edges = np.arange(8,1000)
ells,nls = lensing.get_nl('planck')
fsky = 0.65
lens = pyfisher.get_lensing_fisher(bin_edges,ells,nls,fsky)
F = lens #bao+lens
F.delete(['w0','wa','ok','nnu','tau','mnu'])

deriv_path = "pyfisher/data/v20201120_s8_derivs/"
F = pyfisher.reparameterize(F,['s8','om','H0'],fids,deriv_path)

print(F)

s8 = pyfisher.get_s8(zs=[0.],params=fids)[0]
fids['s8'] = s8
fids['om'] = (fids['omch2'] + fids['ombh2'])/(fids['H0']/100)**2.

print(fids)

pyfisher.contour_plot(F,fids,'contour.png',name=None)
