import numpy as np
import sys,os
import ConfigParser
import argparse
#import LikeRealAL as Lens
#import liteMap
#from astLib import astWCS
from makeDerivs_axionCAMB import getClsAxionCamb


verbose = True
parser = argparse.ArgumentParser(description='Make Derivatives using axionCAMB',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--out_pre",type=str)
parser.add_argument("--spec",type=str)
parser.add_argument("--AccuracyBoost",type=bool)
parser.add_argument("--templateIni",type=str)

parser.add_argument("--paramName",type=str)
parser.add_argument("--stepSize",type=float)
parser.add_argument("--output",type=str)
parser.add_argument("--seed",type=int)

args = parser.parse_args()
out_pre = args.out_pre
spec = args.spec
AccuracyBoost = args.AccuracyBoost
templateIni = args.templateIni
paramName = args.paramName
output = args.output
seed = args.seed

nameMap={'H0':'hubble','YHe':'helium_fraction','nnu':'massless_neutrinos','s_pivot':'pivot_scalar','t_pivot':'pivot_tensor','As':'scalar_amp(1)','ns':'scalar_spectral_index(1)','tau':'re_optical_depth','num_massive_neutrinos':'massive_neutrinos','TCMB':'temp_cmb','lmax':'l_max_scalar','kmax':'k_eta_max_scalar'}
inv_nameMap = {v: k for k, v in nameMap.items()}
if paramName in nameMap:
	paramName = nameMap[paramName]

#Get fiducial
iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/makeDerivs_axionCAMB_input.ini' 
Config = ConfigParser.SafeConfigParser()
Config.optionxform = str
Config.read(iniFile)

fparams = {}
fidscript = ''
for (key, val) in Config.items('axionCAMB'):
        fidscript += key+' = '+val+'\n'
        if key in nameMap:
                key = nameMap[key]
	if ',' in val:
		param, step = val.split(',')
		fparams[key] = float(param)
	else:
		if key in ['l_max_scalar','massive_neutrinos']:
			fparams[key] = int(val)
		else:
			fparams[key] = float(val)
if 'massless_neutrinos' in fparams:
        fparams['massless_neutrinos'] -= fparams['massive_neutrinos']
                    
if not('omnuh2' in fparams) and ('mnu' in fparams):
	fparams['omnuh2'] = round(fparams['mnu']/93.14,6)

# Fid or deriv 
if paramName == '0':
	print "Calculating and saving fiducial cosmology..."
	CLS = getClsAxionCamb(out_pre,spec,templateIni,fparams,AccuracyBoost=AccuracyBoost,seed=seed)
	output = output+'_axion_fCls.csv'
        # Save fid + stepsize
        fidscript ='[camb]\n'+fidscript
        filename = args.output+'_axion_fid.csv'
        with open(filename,'w') as tempFile:
                tempFile.write(fidscript)
                                                
else:
        h = args.stepSize
        print "Calculating forward difference for ", paramName
        pparams = fparams.copy()
        pparams[paramName] = fparams[paramName] + 0.5*h
        if paramName =='mnu':
            pparams['omnuh2'] = round(pparams['mnu']/93.14,6)
        pCls = getClsAxionCamb(out_pre+str(seed),spec,templateIni,pparams,AccuracyBoost=AccuracyBoost,seed=seed)

        print "Calculating backward difference for ", paramName
        mparams = fparams.copy()
        mparams[paramName] = fparams[paramName] - 0.5*h
        if paramName =='mnu':
            mparams['omnuh2'] = round(mparams['mnu']/93.14,6)
        mCls = getClsAxionCamb(out_pre+str(seed),spec,templateIni,mparams,AccuracyBoost=AccuracyBoost,seed=seed)

	CLS = (pCls-mCls)/h
	if paramName in inv_nameMap:
		output = output+'_axion_dCls_'+inv_nameMap[paramName]+'.csv'
	else:
                output = output+'_axion_dCls_'+paramName+'.csv'
np.savetxt(output,CLS,delimiter=",")
