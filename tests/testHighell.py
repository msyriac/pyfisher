from pyfisher.lensInterface import lensNoise
import sys, os
from ConfigParser import SafeConfigParser 

# Read config
iniFile = "input/params.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)


px = 180*60/60000.
gradCut = 60000.

expName = "DM-18arcsec"
lensName = "lensing"
ls,Nls,ellbb,dlbb,efficiency = lensNoise(Config,expName,lensName,px=px,gradCut=gradCut,bigell=gradCut,plot=True)



print efficiency
print ls
print Nls


from orphics.tools.io import Plotter

pl = Plotter(scaleY='log')
pl.add(ls,Nls)
pl._ax.set_ylim(1e-12,1e-6)
#pl._ax.set_xlim(2,4000)
pl.done("nls.png")
