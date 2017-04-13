import glob
import numpy as np
from mmUtils import Plotter
import camb
import sys

# pars = camb.CAMBparams()
# pars.set_cosmology(tau=0.01)
# pars.InitPower.set_params()
# results = camb.get_results(pars)
# powers =results.get_cmb_power_spectra(pars)
# #CL=powers[spec]*1.e12*params['tcmb']**2.
# lensArr = results.get_lens_potential_cls(lmax = 2000)
# clphi = lensArr[:,0]

# pars = camb.CAMBparams()
# pars.set_cosmology(tau=0.09)
# pars.InitPower.set_params()
# results = camb.get_results(pars)
# powers =results.get_cmb_power_spectra(pars)
# #CL=powers[spec]*1.e12*params['tcmb']**2.
# lensArr2 = results.get_lens_potential_cls(lmax = 2000)
# clphi2 = lensArr2[:,0]

# print np.nanmean((clphi-clphi2)*100./clphi)

# sys.exit()



derivRoot = "defAccExt_unlensed_dCls_"


fileList = glob.glob("output/"+derivRoot+"*.csv")
farr = np.loadtxt("output/"+derivRoot[:-5]+"fCls.csv",delimiter=',')

arrs = {}
for fileN in fileList:
    lab = fileN[len("output/"+derivRoot):-4]
    arrs[lab] = np.loadtxt(fileN,delimiter=',')


specList = ['TT','EE','BB','TE','KK','KT']
for spec in specList:
    if spec=='BB': continue
    pls = Plotter(scaleY='log',scaleX='log')

    

    ind = specList.index(spec)
    
    for lab in arrs:
        arr = arrs[lab]
        y = arr[:,ind]**2./farr[:,ind]**2.
        if lab=='tau': ls = '--'
        else: ls = "-"
        pls.add(range(arr.shape[0]),y,label=lab,ls=ls)


    pls.legendOn(loc='upper right',labsize=8)
    pls._ax.set_xlim(20.,4000.)
    pls.done("output/d"+spec+".png")
