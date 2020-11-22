import sys, os
from ConfigParser import SafeConfigParser 
import cPickle as pickle
import numpy as np
from scipy.interpolate import interp1d
import argparse
from pyfisher.lensInterface import lensNoise
from pyfisher.clFisher import tryLoad, calcFisher, loadFishers, noiseFromConfig, rSigma
from orphics.tools.io import dictFromSection, listFromConfig, printC
import orphics.tools.cmb as cmb
from orphics.theory.gaussianCov import LensForecast
import cPickle as pickle
from orphics.tools.io import Plotter

# Get the name of the experiment and lensing type from command line
parser = argparse.ArgumentParser(description='Run a Fisher test.')
parser.add_argument('expName', type=str,help='The name of the experiment in input/params.ini')
parser.add_argument('lensName',type=str,help='The name of the CMB lensing section in input/params.ini. ',default="")
args = parser.parse_args()
expName = args.expName
lensName = args.lensName

# mnu Forecast ============================================

TCMB = 2.7255e6

# Read config
iniFile = "input/params.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)



efficiencies = []
mnus = []
sns = []
rs = []
rdelens = []

if True:
    # Get lensing noise curve. If you want to override something from the Config file in order to make plots varying it,
    # change from None to the value you want.

    ls,Nls,ellbb,dlbb,efficiency = lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None)

    efficiencies.append(efficiency)
    
    printC("Delensing efficiency: "+ str(efficiency) + " %",color="green",bold=True)

    # File root name for Fisher derivatives
    derivRoot = Config.get("fisher","derivRoot")

    # Get list of parameters
    paramList = Config.get("fisher","paramList").split(',')

    # Load fiducials and derivatives
    fidCls = tryLoad(derivRoot+'_fCls.csv',',')
    dCls = {}
    for paramName in paramList:
        dCls[paramName] = tryLoad(derivRoot+'_dCls_'+paramName+'.csv',',')


    # Load other Fisher matrices to add
    otherFisher = loadFishers(Config.get('fisher','otherFishers').split(','))

    # Get CMB noise functions and ell ranges. Note that the same overriding is possible but here the beams and noises have to be lists for the different frequencies.
    fnTT, fnEE = noiseFromConfig(Config,expName,TCMB=TCMB,beamsOverride=None,noisesOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None)
    tellmin,tellmax = listFromConfig(Config,expName,'tellrange')
    pellmin,pellmax = listFromConfig(Config,expName,'pellrange')

    # Pad CMB lensing noise with infinity outside L ranges
    kellmin,kellmax = listFromConfig(Config,'lensing','Lrange')
    fnKK = cmb.noise_pad_infinity(interp1d(ls,Nls,fill_value=np.inf,bounds_error=False),kellmin,kellmax)

    # Decide on what ell range to calculate the Fisher matrix
    ellrange = np.arange(min(tellmin,pellmin,kellmin),max(tellmax,pellmax,kellmax)).astype(int)
    # Get fsky
    fsky = Config.getfloat(expName,'fsky')
    # Calculate the Fisher matrix and add to other Fishers
    # Fisher = otherFisher+calcFisher(paramList,ellrange,fidCls,dCls,fnTT,fnEE,fnKK,fsky,verbose=True)

    # # Get prior sigmas and add to Fisher
    # priorList = Config.get("fisher","priorList").split(',')
    # for prior,param in zip(priorList,paramList):
    #     try:
    #         priorSigma = float(prior)
    #     except ValueError:
    #         continue
    #     ind = paramList.index(param)
    #     Fisher[ind,ind] += 1./priorSigma**2.

    # # get index of mnu and print marginalized constraint
    # indMnu = paramList.index('mnu')
    # mnu = np.sqrt(np.linalg.inv(Fisher)[indMnu,indMnu])*1000.
    # printC("Sum of neutrino masses 1-sigma: "+ str(mnu) + " meV",color="green",bold=True)

    # mnus.append(mnu)
    # CLKK S/N ============================================

    # Calculate Clkk S/N
    Clkk = fidCls[:,4]
    frange = np.array(range(len(Clkk)))
    snrange = np.arange(kellmin,kellmax)
    LF = LensForecast()
    LF.loadKK(frange,Clkk,ls,Nls)
    sn,errs = LF.sn(snrange,fsky,"kk")
    printC("Lensing autopower S/N: "+ str(sn),color="green",bold=True)


    sns.append(sn)
    




