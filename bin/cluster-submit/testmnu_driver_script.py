import camb
from camb import model, initialpower
import numpy as np
import sys
import ConfigParser
import os
import argparse
from driver import FisherForecast

verbose = True

parser = argparse.ArgumentParser(description='Test Fiducial for Fisher Code - driver',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--testParam",type=str)
parser.add_argument("--output",type=str)
args = parser.parse_args()

testParam = args.testParam
output = args.output

# Read config
IniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testmnu_driver_input.ini'
Config = ConfigParser.SafeConfigParser()
Config.optionxform = str
Config.read(IniFile)
paramList = Config.get('general','paramList').split(',')

iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testmnu_input_optimal.ini'
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

testRange = {}
for (key, val) in config.items('testmnu'):
        pmin,pmax,numstep = [float(x) for x in val.split(',')]
        pcenter = fparams[key]
        mstep = (pcenter-pmin)/(numstep/2)
        pstep = (pmax-pcenter)/(numstep/2)
        testRange[key] = np.append(np.arange(pmin,pcenter,mstep),np.arange(pcenter,pmax+pstep/2,pstep))
        # round the value to about 5 significant digits                                                                                                      
        decimal = int(-(np.ceil(np.log10(pmin))-5))
        testRange[key] = [round(x,decimal) for x in testRange[key]]

print testParam,testRange[testParam]
# FISHER
data = np.zeros([len(testRange[testParam]),(1+len(paramList))])
i = 0
for testVal in testRange[testParam]:
        prefix = testParam+'_'+str(testVal)
        F = FisherForecast(IniFile,prefix=prefix)
        F.calcFisher(verbose = True)
        j = 0
        data[i,j]= testVal
        j += 1
        for paramName in paramList:
                data[i,j] = F.margSigma(paramName)
                j += 1
        i+=1
np.savetxt(output+'_'+testParam+'.csv',data,delimiter=",")

#np.savetxt(os.environ['FISHER_DIR']+"/output/June3_testfid_mnu_"+testParam+".csv",mnudata,delimiter=",")
#np.savetxt(os.environ['FISHER_DIR']+"/output/June3_testfid_H0_"+testParam+".csv",H0data,delimiter=",")
