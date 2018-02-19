import sys, os
from configparser import SafeConfigParser 
import cPickle as pickle
import numpy as np
from scipy.interpolate import interp1d
import argparse
from pyfisher.lensInterface import lensNoise
from pyfisher.clFisher import tryLoad, calcFisher, loadFishers, noiseFromConfig, rSigma
from orphics.io import dict_from_section, list_from_config, cprint, Plotter
from orphics.cosmology import LensForecast
import cPickle as pickle
from orphics import cosmology

# Get the name of the experiment and lensing type from command line
parser = argparse.ArgumentParser(description='Run a Fisher test.')
parser.add_argument('expName', type=str,help='The name of the experiment in input/params.ini')
parser.add_argument('lensName',type=str,help='The name of the CMB lensing section in input/params.ini. ')
parser.add_argument('saveName',type=str,help='Name of file to save Fisher to.')
args = parser.parse_args()
expName = args.expName
lensName = args.lensName
saveName = args.saveName


TCMB = 2.7255e6

# Read config
iniFile = "input/params.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)

ls,Nls,ellbb,dlbb,efficiency,cc = lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None)

# File root name for Fisher derivatives
derivRoot = Config.get("fisher","derivRoot")

# Get list of parameters
paramList = Config.get("fisher","paramList").split(',')

# Load fiducials and derivatives
fidCls = tryLoad(derivRoot+'_fCls.csv',',')
dCls = {}
for paramName in paramList:
    dCls[paramName] = tryLoad(derivRoot+'_dCls_'+paramName+'.csv',',')


# Get CMB noise functions and ell ranges. 
fnTT, fnEE = noiseFromConfig(Config,expName,TCMB=TCMB,beamsOverride=None,noisesOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None)

tellmin,tellmax = list_from_config(Config,expName,'tellrange')
pellmin,pellmax = list_from_config(Config,expName,'pellrange')

# Pad CMB lensing noise with infinity outside L ranges
kellmin,kellmax = list_from_config(Config,'lensing','Lrange')
fnKK = cosmology.noise_pad_infinity(interp1d(ls,Nls,fill_value=np.inf,bounds_error=False),kellmin,kellmax)

# Decide on what ell range to calculate the Fisher matrix
ellrange = np.arange(min(tellmin,pellmin,kellmin),max(tellmax,pellmax,kellmax)).astype(int)
# Get fsky
fsky = Config.getfloat(expName,'fsky')
# Calculate the Fisher matrix and add to other Fishers
Fisher = calcFisher(paramList,ellrange,fidCls,dCls,lambda x: fnTT(x)*TCMB**2.,lambda x: fnEE(x)*TCMB**2.,fnKK,fsky,verbose=True)

np.savetxt(saveName,Fisher)
