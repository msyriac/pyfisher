import sys, os
from ConfigParser import SafeConfigParser 
import numpy as np
from scipy.interpolate import interp1d
import argparse
from pyfisher.clFisher import noiseFromConfig
from orphics.io import list_from_config
import orphics.cosmology as cosmo
from orphics.io import Plotter

# Get the name of the experiment and lensing type from command line
parser = argparse.ArgumentParser(description='Run a Fisher test.')
parser.add_argument('expName', type=str,help='The name of the experiment in input/params.ini')
parser.add_argument('lensName',type=str,help='The name of the CMB lensing section in input/params.ini. ',default="")
parser.add_argument('saveName', nargs='?',type=str,help='Suffix for plots ',default="")
args = parser.parse_args()
expName = args.expName
lensName = args.lensName
saveName = args.saveName

# mnu Forecast ============================================

TCMB = 2.7255e6

# Read config
iniFile = "input/params_testFixedTime_SO.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)

# Fixed time
fskyList = np.hstack([0.02,np.arange(0.05,0.55,0.05)])
noiseList = 2.0*np.sqrt(fskyList/0.1)

for noiseNow,fskyNow in zip(noiseList,fskyList):
    
    #print np.round(noiseNow,2),fskyNow
    
    # Get CMB noise functions and ell ranges. Note that the same overriding is possible but here the beams and noises have to be lists for the different frequencies.
    fnTT, fnEE = noiseFromConfig(Config,expName,TCMB=TCMB,beamsOverride=None,noisesOverride=[noiseNow],lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None,tellminOverride=0.,pellminOverride=0.)
    
    tellmin,tellmax = list_from_config(Config,expName,'tellrange')
    pellmin,pellmax = list_from_config(Config,expName,'pellrange')

    nTT = np.vstack([np.arange(tellmax),fnTT(np.arange(tellmax))]).T
    nEE = np.vstack([np.arange(pellmax),fnEE(np.arange(pellmax))]).T

    np.savetxt('output/'+saveName+'_whiteNoise_noatm_TT_noiseT_'+str(np.round(noiseNow,2))+'.csv',nTT)
    np.savetxt('output/'+saveName+'_whiteNoise_noatm_EE_noiseT_'+str(np.round(noiseNow,2))+'.csv',nEE)
