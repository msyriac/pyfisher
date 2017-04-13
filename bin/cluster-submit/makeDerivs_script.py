import camb
from camb import model, initialpower
import numpy as np
import sys
import ConfigParser
import os
import argparse
#import LikeRealAL as Lens
#import liteMap
#from astLib import astWCS
from makeDerivs import getHubbleCosmology,getPowerCamb
verbose = True

parser = argparse.ArgumentParser(description='Make Derivatives',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--spec",type=str)
parser.add_argument("--out_pre",type=str)
parser.add_argument("--AccuracyBoost",type=bool)


parser.add_argument("--paramName",type=str)
parser.add_argument("--stepSize",type=float)
parser.add_argument("--output",type=str)
args = parser.parse_args()
spec = args.spec
AccuracyBoost = args.AccuracyBoost
paramName = args.paramName
stepSize = args.stepSize

#Get fiducial
iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/makeDerivs_input.ini' 
Config = ConfigParser.SafeConfigParser()
Config.optionxform = str
Config.read(iniFile)

fparams = {}
fidscript = ''
for (key, val) in Config.items('camb'):
        if ',' in val:
                param, step = val.split(',')
                fparams[key] = float(param)
        else:
                fparams[key] = float(val)
	fidscript += key+' = '+val+'\n'
if not('H0' in fparams):
	fparams['H0'] = getHubbleCosmology(fparams['theta100'],fparams)
try:
        derivForm = Config.getint('general','derivForm')
except:
        derivForm = 0
            
# Fid or deriv 
if paramName == '0':
        CLS = getPowerCamb(fparams,spec,AccuracyBoost=AccuracyBoost)
	output = args.output+'_fCls.csv'
	# Save fid + stepsize
	fidscript ='[camb]\n'+fidscript
	filename = args.output+'_fid.csv'
	with open(filename,'w') as tempFile:
		tempFile.write(fidscript)

else:
        if derivForm == 0:
                h = stepSize
                pparams = fparams.copy()
                pparams[paramName] = fparams[paramName] + 0.5*h
                if paramName=='theta100':
                        pparams['H0'] = getHubbleCosmology(theta=pparams['theta100'],params=pparams)
                pCls = getPowerCamb(pparams,spec,AccuracyBoost=AccuracyBoost)
                
                mparams = fparams.copy()
                mparams[paramName] = fparams[paramName] - 0.5*h
                if paramName=='theta100':
                        mparams['H0'] = getHubbleCosmology(theta=mparams['theta100'],params=mparams)
                mCls = getPowerCamb(mparams,spec,AccuracyBoost=AccuracyBoost)
                CLS = (pCls-mCls)/h
        elif derivForm == 1:
                h = 0.5*stepSize
                
                params1 = fparams.copy()
                params2 = fparams.copy()
                params3 = fparams.copy()
                params4 = fparams.copy()
                params1[paramName] = fparams[paramName] + 2.*h
                params2[paramName] = fparams[paramName] + h
                params3[paramName] = fparams[paramName] - h
                params4[paramName] = fparams[paramName] - 2.*h
                if paramName=='theta100':
                        params1['H0'] = getHubbleCosmology(theta=params1['theta100'],params=params1)
                        params2['H0'] = getHubbleCosmology(theta=params2['theta100'],params=params2)
                        params3['H0'] = getHubbleCosmology(theta=params3['theta100'],params=params3)
                        params4['H0'] = getHubbleCosmology(theta=params4['theta100'],params=params4)
                Cls1 = getPowerCamb(params1,spec,AccuracyBoost=AccuracyBoost)
                Cls2 = getPowerCamb(params2,spec,AccuracyBoost=AccuracyBoost)
                Cls3 = getPowerCamb(params3,spec,AccuracyBoost=AccuracyBoost)
                Cls4 = getPowerCamb(params4,spec,AccuracyBoost=AccuracyBoost)
                
                dCls = (- Cls1 + 8.*Cls2 - 8.*Cls3 + Cls4)/(12.*h)
                
	output = args.output+'_dCls_'+paramName+'.csv'
np.savetxt(output,CLS,delimiter=",")
