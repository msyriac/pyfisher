import camb
from camb import model, initialpower
import numpy as np
import sys
import ConfigParser
import os
import argparse
from driver import FisherForecast

verbose = True

parser = argparse.ArgumentParser(description='Test Step Size for Fisher Code - driver',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--testParam",type=str)
parser.add_argument("--testVal",type=float)
parser.add_argument("--output",type=str)
parser.add_argument("--index",type=int)
args = parser.parse_args()

testParam = args.testParam
testVal = args.testVal
output = args.output
index = args.index

# Read config
IniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testStepSize_driver_input.ini'
Config = ConfigParser.SafeConfigParser()
Config.optionxform = str
Config.read(IniFile)
paramList = Config.get('general','paramList').split(',')

iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testStepSize_input.ini'
config = ConfigParser.SafeConfigParser()
config.optionxform = str
config.read(iniFile)

fparams = {}
for (key, val) in config.items('camb'):
        if ',' in val:
                param, step= val.split(',')
                fparams[key] = float(param)
        else:
                fparams[key] = float(val)
'''
testRange = {}
logspace = config.getboolean('testStepSize','logspace')
pmin,pmax,numstep = [float(x) for x in config.get('testStepSize','testRange').split(',')]
if logspace:
        fraction = np.logspace(np.log10(pmin),np.log10(pmax),numstep)
	order = 2
else:
        fraction = np.linspace(pmin,pmax,numstep)
	order = 3

testRange[testParam]=abs(fraction*fparams[testParam])
# round the value to about 'order' significant digits                                                                                                          
decimal = int(-(np.ceil(np.log10(testRange[testParam][0]))-order))
testRange[testParam] = [round(x,decimal) for x in testRange[testParam]]
print testRange
'''
# FISHER
data = np.zeros([1,(1+len(paramList))])

dClsRoot={}
dClsRoot[testParam] = output+'_'+testParam+'_'+str(testVal)
F = FisherForecast(IniFile,dClsRoot=dClsRoot)
F.calcFisher(verbose=True)
j = 0
data[0,j]= testVal
j += 1
for paramName in paramList:
        data[0,j] = F.margSigma(paramName)
        j += 1
postfix = '_LCDM+mnu_S4_'
np.savetxt(output+'_'+testParam+postfix+'part'+str(index)+'.csv',data,delimiter=",")
