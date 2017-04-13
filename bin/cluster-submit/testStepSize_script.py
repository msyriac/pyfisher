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

parser = argparse.ArgumentParser(description='Test Step Size for Fisher Code',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--spec",type=str)
parser.add_argument("--out_pre",type=str)
parser.add_argument("--AccuracyBoost",type=bool)

parser.add_argument("--testParam",type=str)
parser.add_argument("--testVal",type=float)
parser.add_argument("--output",type=str)
args = parser.parse_args()

#Get fiducial
iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testStepSize_input.ini'
config = ConfigParser.SafeConfigParser()
config.optionxform = str
config.read(iniFile)

fparams = {}
for (key, val) in config.items('camb'):
        if ',' in val:
                param, step = val.split(',')
                fparams[key] = float(param)
        else:
                fparams[key] = float(val)
if not('H0' in fparams):
        fparams['H0'] = getHubbleCosmology(fparams['theta100'],fparams)
# Assigning testVal(stepsize) to testParam
spec = args.spec
AccuracyBoost = args.AccuracyBoost
paramName = args.testParam
h = args.testVal

print "Calculating forward difference for ", paramName
pparams = fparams.copy()
pparams[paramName] = fparams[paramName] + 0.5*h
if paramName=='theta100':
	pparams['H0'] = getHubbleCosmology(theta=pparams['theta100'],params=pparams)
pCls = getPowerCamb(pparams,spec,AccuracyBoost=AccuracyBoost)

print "Calculating backward difference for ", paramName
mparams = fparams.copy()
mparams[paramName] = fparams[paramName] - 0.5*h
if paramName=='theta100':
	mparams['H0'] = getHubbleCosmology(theta=mparams['theta100'],params=mparams)
mCls = getPowerCamb(mparams,spec,AccuracyBoost=AccuracyBoost)

dCls = (pCls-mCls)/h

np.savetxt(args.output,dCls,delimiter=",")
