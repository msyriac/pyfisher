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
parser.add_argument("-t", "--tt",     type=str,  default=None,help="Dimensionless beam-deconvolved TT noise curve file to override experiment.")
parser.add_argument("-p", "--pp",     type=str,  default=None,help="Dimensionless beam-deconvolved PP (EE/BB) noise curve file to override experiment.")
args = parser.parse_args()
expName = args.expName
lensName = args.lensName
saveName = args.saveName


try:
    elltt,ntt = np.loadtxt(args.tt,usecols=[0,1],unpack=True)
    noise_func_tt = interp1d(elltt,ntt,bounds_error=False,fill_value=np.inf)
except:
    noise_func_tt = None

try:
    ellee,nee = np.loadtxt(args.pp,usecols=[0,1],unpack=True)
    noise_func_ee = interp1d(ellee,nee,bounds_error=False,fill_value=np.inf)
except:
    noise_func_ee = None
    
    
TCMB = 2.7255e6

# Read config
iniFile = "input/params_local.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)

ls,Nls,ellbb,dlbb,efficiency,cc = lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None,noiseFuncT=noise_func_tt,noiseFuncP=noise_func_ee)

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
if (noise_func_tt is None) or (noise_func_ee is None):
    fnTT, fnEE = noiseFromConfig(Config,expName,TCMB=TCMB,beamsOverride=None,noisesOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None)

tellmin,tellmax = list_from_config(Config,expName,'tellrange')
pellmin,pellmax = list_from_config(Config,expName,'pellrange')

if (noise_func_tt is not None):
    fnTT = cosmology.noise_pad_infinity(noise_func_tt,tellmin,tellmax)
if (noise_func_ee is not None):
    fnEE = cosmology.noise_pad_infinity(noise_func_ee,pellmin,pellmax)
    

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
