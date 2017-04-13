import sys
import numpy as np
from driver import FisherForecast
import ConfigParser
import matplotlib.pyplot as plt
import os
from gpcInterface import jobMaker

#iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testmnu_input.ini'
iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testmnu_input_optimal.ini'
config = ConfigParser.SafeConfigParser()
config.optionxform = str
config.read(iniFile)

numCoresAsked = config.getint('general','numCores')
walltime = config.get('general','walltime')
saveRoot = config.get('general', 'saveRoot')
spec = config.get('general','spec')
out_pre = config.get('general','output_prefix')
AccuracyBoost = config.getboolean('general','AccuracyBoost')
queue = config.get('general','queue')
projectName = config.get('general','projectName')

fparams = {}
for (key, val) in config.items('camb'):
        if ',' in val:
                param, step= val.split(',')
                fparams[key] = float(param)
        else:
                fparams[key] = float(val)

#testList = config.get('testStepSize','testList').split(',')
testList = []
for (key, val) in config.items('testmnu'):
        testList.append(key)

# a magical formula                                                                                                                                            
N = 0
N = len(testList)
#N = 1
print "Number of tasks: ", N
tasksPerCore = np.ceil(float(N)/float(numCoresAsked))
numCores = float(N) / tasksPerCore
print "Tasks per core: ",tasksPerCore
print "Needed cores: ",int(np.ceil(numCores)),"/",numCoresAsked
print "Needed nodes: ",int(np.ceil(numCores/8)),"/",numCoresAsked/8
#sys.exit()

commandPreFix = 'source ~/.bashrc'
Jobs = jobMaker(projectName=projectName,commandPreFix=commandPreFix,numCores=8,queue=queue,walltime=walltime,jobRoot='/gpfs/scratch/nhnguyen/testmnuDump/jobs/')
#optionsCommon = ' --spec ' + spec + ' --out_pre ' + out_pre + ' --AccuracyBoost ' + str(AccuracyBoost)
optionsCommon = ''
i = 0
nowCommand = ''
#prefix = ''
for testParam in testList:
        #prefix = testParam+'_'+str(testVal)
        i+=1
        saveSuffix = out_pre+spec
        options = optionsCommon + ' --testParam ' + testParam + ' --output ' + saveRoot+saveSuffix
        nowCommand += ' python '+os.environ['FISHER_DIR']+'/bin/cluster-submit/testmnu_driver_script.py '+options+';'
        if i%tasksPerCore == 0 or i==N:
                Jobs.addJob(nowCommand)
                nowCommand = ''
#print len(Jobs.scripts),Jobs.scripts 
Jobs.submit()
