import os, sys
import glob
import time
import subprocess
import ConfigParser
import numpy as np
from gpcInterface import jobMaker
    
iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/oldtestStepSize_input.ini'
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

#paramList = []
fparams = {}
for (key, val) in config.items('camb'):
        #paramList.append(key)
        fparams[key]=float(val)

testRange = {}
testList = []

for (key, val) in config.items('oldtestStepSize'):
        pmin,pmax,numstep = [float(x) for x in val.split(',')]
	testList.append(key)
        testRange[key] = np.linspace(pmin,pmax,numstep)
	# round the value to about 6 significant digits
	decimal = int(-(np.ceil(np.log10(abs(pmin)))-6))
        testRange[key] = [round(x,decimal) for x in testRange[key]]
	#print key,pmin,pcenter,pmax

# a magical formula
numCores=2
N = 0
for testParam in testList:
        # for each value in testRange, needs 1 call for fid
        N+=len(testRange[testParam])
	print testParam,testRange[testParam]
print "Number of tasks: ", N
print "Number of cores running per node:",numCores,"/",8
tasksPerCore = np.ceil(float(N)/float(numCoresAsked)*8/numCores)
numCoresNeed = float(N) / tasksPerCore
print "Tasks per core: ",int(tasksPerCore)
print "Needed cores: ",int(np.ceil(numCoresNeed)),"/",numCoresAsked
print "Needed nodes: ",int(np.ceil(numCoresNeed/numCores)),"/",numCoresAsked/8
#sys.exit()

commandPreFix = 'source ~/.bashrc'
Jobs = jobMaker(projectName=projectName,commandPreFix=commandPreFix,numCores=numCores,queue=queue,walltime=walltime,jobRoot='/gpfs/scratch/nhnguyen/oldtestStepSize/jobs/')
optionsCommon = ' --spec ' + spec + ' --out_pre ' + out_pre + ' --AccuracyBoost ' + str(AccuracyBoost)
i = 0
nowCommand = ''
prefix = ''
for testParam in testList:
        for testVal in testRange[testParam]:
                prefix = testParam+'_'+str(testVal)
                # Save fiducials
                i+=1
                saveSuffix = out_pre+spec+"_"+prefix+"_fCls.csv"
                options = optionsCommon + ' --testParam ' + testParam + ' --testVal ' + str(testVal) + ' --output ' + saveRoot+saveSuffix 
                nowCommand += ' python '+os.environ['FISHER_DIR']+'bin/cluster-submit/oldtestStepSize_script.py '+options+';'
		if i%tasksPerCore == 0 or i==N:
			Jobs.addJob(nowCommand)
                        nowCommand = ''
#print len(Jobs.scripts),Jobs.scripts
Jobs.submit()
