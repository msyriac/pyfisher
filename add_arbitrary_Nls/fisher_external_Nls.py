import sys, os
from ConfigParser import SafeConfigParser 
import cPickle as pickle
import numpy as np
from scipy.interpolate import interp1d
import argparse
from pyfisher.lensInterface import lensNoise
from pyfisher.clFisher import tryLoad, calcFisher, loadFishers, noiseFromConfig, rSigma
from orphics.io import list_from_config, cprint
import orphics.cosmology as cosmo
#from orphics.cosmology import LensForecast
from orphics.io import Plotter

# Get the name of the experiment and lensing type from command line
parser = argparse.ArgumentParser(description='Run a Fisher test.')
parser.add_argument('expName', type=str,help='The name of the experiment in input/params.ini')
parser.add_argument('lensName',type=str,help='The name of the CMB lensing section in input/params.ini. ',default="")
#parser.add_argument('saveName', nargs='?',type=str,help='Suffix for plots ',default="")
parser.add_argument('saveName',type=str,help='Suffix for plots ',default="")
parser.add_argument('noiseFilePrefix',nargs='?', type=str,help='Prefix of external noise files',default=None)

args = parser.parse_args()
expName = args.expName
lensName = args.lensName
saveName = args.saveName
noiseFile = args.noiseFilePrefix
print expName,lensName,saveName,noiseFile
# mnu Forecast ============================================

TCMB = 2.7255e6

# Read config
iniFile = "input/params_testFixedTime_SO.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)

fskyList = np.hstack([0.02,np.arange(0.05,0.55,0.05)])
noiseList = 2.0*np.sqrt(fskyList/0.1)

efficiencies = []
mnus = []
sns = []


def getNoise
for noiseNow,fskyNow in zip(noiseList,fskyList):

    # Get CMB noise functions from external files
    if noiseFile is not None:
        ellT, nlTT = np.loadtxt(noiseFile+'_TT_noiseT_'+str(np.round(noiseNow,2))+'.csv',unpack=True)
        fnTT = interp1d(ellT,nlTT,bounds_error=False,fill_value=np.inf)
        ellP, nlEE = np.loadtxt(noiseFile+'_EE_noiseT_'+str(np.round(noiseNow,2))+'.csv',unpack=True)
        fnEE = interp1d(ellP,nlEE,bounds_error=False,fill_value=np.inf)
        funcNoise = True
        print 'Using noise files '+noiseFile+'_TT_noiseT_'+str(np.round(noiseNow,2))+'.csv'
    else:
        fnTT, fnEE = noiseFromConfig(Config,expName,TCMB=TCMB,beamsOverride=None,noisesOverride=[noiseNow],lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None)
        funcNoise = False
        print 'Not using noise files'

    # Get lensing noise curve. If you want to override something from the Config file in order to make plots varying it,
    # change from None to the value you want.
    if funcNoise:
        ls,Nls,ellbb,dlbb,efficiency,cc = lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=noiseNow,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None,noiseFuncT=lambda x: fnTT(x),noiseFuncP=lambda x: fnEE(x))
    else:
        ls,Nls,ellbb,dlbb,efficiency,cc = lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=noiseNow,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None)
        
    efficiencies.append(efficiency)
    
    cprint("Delensing efficiency: "+ str(efficiency) + " %",color="green",bold=True)

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
    try:
        otherFisher = loadFishers(Config.get('fisher','otherFishers').split(','))
    except:
        otherFisher = 0.

    tellmin,tellmax = list_from_config(Config,expName,'tellrange')
    pellmin,pellmax = list_from_config(Config,expName,'pellrange')

    '''
    # Uncomment this when checking with white noise
    # Reason: 
    # - white noise files need to go to l=16000 since in modLMap, lmax~15500 
    # - when fnTT/EE is fed to Fisher, re-establish the proper bounds from tellmin,tellmax

    tell = np.arange(tellmin,tellmax)
    pell = np.arange(pellmin,pellmax)

    fnTT = interp1d(tell,fnTT(tell),bounds_error=False,fill_value=np.inf)
    fnEE = interp1d(pell,fnEE(pell),bounds_error=False,fill_value=np.inf)
    '''
    # Pad CMB lensing noise with infinity outside L ranges
    kellmin,kellmax = list_from_config(Config,'lensing','Lrange')
    fnKK = cosmo.noise_pad_infinity(interp1d(ls,Nls,fill_value=np.inf,bounds_error=False),kellmin,kellmax)

    # Decide on what ell range to calculate the Fisher matrix
    ellrange = np.arange(min(tellmin,pellmin,kellmin),max(tellmax,pellmax,kellmax)).astype(int)
    # Get fsky
    fsky = fskyNow #Config.getfloat(expName,'fsky')
    # Calculate the Fisher matrix and add to other Fishers
    Fisher = otherFisher+calcFisher(paramList,ellrange,fidCls,dCls,lambda x: fnTT(x)*TCMB**2.,lambda x: fnEE(x)*TCMB**2.,fnKK,fsky,verbose=True)

    # Get prior sigmas and add to Fisher
    priorList = Config.get("fisher","priorList").split(',')
    for prior,param in zip(priorList,paramList):
        try:
            priorSigma = float(prior)
        except ValueError:
            continue
        ind = paramList.index(param)
        Fisher[ind,ind] += 1./priorSigma**2.

    # get index of mnu and print marginalized constraint
    indMnu = paramList.index('mnu')
    mnu = np.sqrt(np.linalg.inv(Fisher)[indMnu,indMnu])*1000.
    cprint("Sum of neutrino masses 1-sigma: "+ str(mnu) + " meV",color="green",bold=True)

    mnus.append(mnu)
    # CLKK S/N ============================================

    # Calculate Clkk S/N
    Clkk = fidCls[:,4]
    frange = np.array(range(len(Clkk)))
    snrange = np.arange(kellmin,kellmax)
    LF = cosmo.LensForecast()
    LF.loadKK(frange,Clkk,ls,Nls)
    sn,errs = LF.sn(snrange,fsky,"kk")
    cprint("Lensing autopower S/N: "+ str(sn),color="green",bold=True)

    sns.append(sn)

outDir = 'output/'+saveName+"_"

np.savetxt(outDir+'mnu.csv',np.vstack([fskyList,mnus]).T)
np.savetxt(outDir+'sn.csv',np.vstack([fskyList,sns]).T)
