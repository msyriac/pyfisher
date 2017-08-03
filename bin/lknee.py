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
parser.add_argument('TorP', type=str,help='T for varying lknee_T, P for varying lknee_P')
parser.add_argument('saveName', nargs='?',type=str,help='Suffix for plots ',default="")
args = parser.parse_args()
expName = args.expName
lensName = args.lensName
saveName = args.saveName
TorP = args.TorP

# mnu Forecast ============================================

TCMB = 2.7255e6

# Read config
iniFile = "input/params.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)


#fskyList = np.append(np.array([0.0025]),np.arange(0.05,0.45,0.05))
fskyList = np.logspace(np.log10(0.001),np.log10(0.7),10)
noiseList = 5.2*np.sqrt(fskyList/0.4)


if TorP=="T":

    outDir = saveName+"_T_"
    lknee_list = [0,1000,2000,3000,4000,5000]

elif TorP=="P":

    outDir = saveName+"_P_"
    lknee_list = [0,200,400,600,800,1000]

else:
    raise ValueError

pl = Plotter(labelX="$f_{\\mathrm{sky}}$",labelY="Clkk S/N")

for lknee in lknee_list:
    efficiencies = []
    mnus = []
    sns = []
    rs = []
    rdelens = []

    for noiseNow,fskyNow in zip(noiseList,fskyList):
        # Get lensing noise curve. If you want to override something from the Config file in order to make plots varying it,
        # change from None to the value you want.

        if TorP=="T":

            ls,Nls,ellbb,dlbb,efficiency,cc = lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=noiseNow,lkneeTOverride=lknee,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None)

        elif TorP=="P":
            
            ls,Nls,ellbb,dlbb,efficiency,cc = lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=noiseNow,lkneeTOverride=None,lkneePOverride=lknee,alphaTOverride=None,alphaPOverride=None)



        efficiencies.append(efficiency)

        printC("Delensing efficiency: "+ str(efficiency) + " %",color="green",bold=True)

        # File root name for Fisher derivatives

        # CLKK S/N ============================================

        # Calculate Clkk S/N
        #Clkk = fidCls[:,4]
        kellmin,kellmax = listFromConfig(Config,'lensing','Lrange')
        fsky = fskyNow #Config.getfloat(expName,'fsky')

        frange = np.arange(0,kellmax) #np.array(range(len(Clkk)))
        Clkk = cc.theory.gCl("kk",frange)
        snrange = np.arange(kellmin,kellmax)
        LF = LensForecast()
        LF.loadKK(frange,Clkk,ls,Nls)
        sn,errs = LF.sn(snrange,fsky,"kk")
        printC("Lensing autopower S/N: "+ str(sn),color="green",bold=True)

        # pl = Plotter(scaleY='log',scaleX='log')
        # pl.add(frange,Clkk)
        # pl.add(ls,Nls)
        # pl._ax.set_ylim(-max(Clkk),max(Clkk))
        # pl.done("clkk.png")



        sns.append(sn)


    pl.add(fskyList,sns,label="lknee_"+TorP+"="+str(lknee))
pl.legendOn(loc="lower right",labsize=10)
pl.done(outDir + "sns.png")


