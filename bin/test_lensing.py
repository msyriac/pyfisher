from __future__ import print_function
from orphics import maps,io,cosmology,stats,lensing
from pixell import enmap
import numpy as np
import os,sys
import pyfisher

fids = pyfisher.get_fiducials()
bao = pyfisher.get_saved_fisher('desi_bao')
lcmb = pyfisher.get_saved_fisher('planck_lowell',0.65)
hcmb = pyfisher.get_saved_fisher('planck_highell',0.65)

bin_edges = np.arange(8,2000)
ells,nls = lensing.get_nl('s4')
fsky = 0.4
lens = pyfisher.get_lensing_fisher(bin_edges,ells,nls,fsky)
F = lcmb+hcmb+bao+lens
F.delete(['w0','wa','ok','nnu'])
F.add_prior('tau',0.01)
print(F.sigmas())

pyfisher.contour_plot(F,fids,'contour.png',name=None)
