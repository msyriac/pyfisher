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

parser = argparse.ArgumentParser(description='Old Test Step Size for Fisher Code',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--spec",type=str)
parser.add_argument("--out_pre",type=str)
parser.add_argument("--AccuracyBoost",type=bool)


parser.add_argument("--testParam",type=str)
parser.add_argument("--testVal",type=float)
parser.add_argument("--output",type=str)
args = parser.parse_args()

#Get fiducial
iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/oldtestStepSize_input.ini' 
Config = ConfigParser.SafeConfigParser()
Config.optionxform = str
Config.read(iniFile)

fparams = {}
for (key, val) in Config.items('camb'):
        if ',' in val:
                param, step = val.split(',')
                fparams[key] = float(param)
        else:
                fparams[key] = float(val)
if not('H0' in fparams):
	fparams['H0'] = getHubbleCosmology(fparams['theta100'],fparams)
# Assigning testVal to testParam
fparams[args.testParam]=args.testVal
if args.testParam=='theta100':
	fparams['H0'] = getHubbleCosmology(fparams['theta100'],fparams)
print "Testing ",args.testParam,"=",args.testVal
CLS = getPowerCamb(fparams,args.spec,AccuracyBoost=args.AccuracyBoost)
np.savetxt(args.output,CLS,delimiter=",")
