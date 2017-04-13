import os, sys
import glob
import time
import subprocess
import ConfigParser
import numpy as np
from gpcInterface import jobMaker
    
iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testStepSize_input.ini'
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
#stepSizes = {}
for (key, val) in config.items('camb'):
	if ',' in val:
		param, step= val.split(',')
		fparams[key] = float(param)
		#stepSizes[key] = float(step)
	else:
		fparams[key] = float(val) 
testList = []
testRange = {}

testList = config.get('testStepSize','testList').split(',')
logspace = config.getboolean('testStepSize','logspace')	
pmin,pmax,numstep = [float(x) for x in config.get('testStepSize','testRange').split(',')]
if logspace:
	fraction = np.logspace(np.log10(pmin),np.log10(pmax),numstep)
	order = 2
else:
	fraction = np.linspace(pmin,pmax,numstep)
	order = 3
for testParam in testList:
	testRange[testParam]=abs(fraction*fparams[testParam])
	# round the value to about 'order' significant digits                                            
	decimal = int(-(np.ceil(np.log10(testRange[testParam][0]))-order))
	testRange[testParam] = [round(x,decimal) for x in testRange[testParam]]

# a magical formula
numCores=1
N = 0
for testParam in testList:
	N += len(testRange[testParam])
	print testParam,testRange[testParam]
print "Number of tasks: ", N
print "Number of cores running per node:",numCores,"/",8
tasksPerCore = np.ceil(float(N)/float(numCoresAsked)*8/numCores)
numCoresNeed = float(N) / tasksPerCore
print "Tasks per core: ",int(tasksPerCore)
print "Needed cores: ",int(np.ceil(numCoresNeed)),"/",numCoresAsked
print "Needed nodes: ",int(np.ceil(numCoresNeed/numCores)),"/",numCoresAsked/8
sys.exit()                               

commandPreFix = 'source ~/.bashrc'
Jobs = jobMaker(projectName=projectName,commandPreFix=commandPreFix,numCores=numCores,queue=queue,walltime=walltime,jobRoot='/gpfs/scratch/nhnguyen/testStepSizeDump/jobs/')
optionsCommon = ' --spec ' + spec + ' --out_pre ' + out_pre + ' --AccuracyBoost ' + str(AccuracyBoost)
i = 0
nowCommand = ''
prefix = ''
for testParam in testList:
	for testVal in testRange[testParam]:
		prefix = testParam+'_'+str(testVal)
		i+=1
		saveSuffix = out_pre+spec+"_"+prefix+"_dCls_"+testParam+".csv"
		options = optionsCommon + ' --testParam ' + testParam + ' --testVal ' + str(testVal) + ' --output ' + saveRoot+saveSuffix 
		nowCommand += ' python '+os.environ['FISHER_DIR']+'/bin/cluster-submit/testStepSize_script.py '+options+';'
		#print nowCommand
		#print i%tasksPerCore, " and ",i," and ", tasksPerCore
		if i%tasksPerCore == 0 or i==N:
                        Jobs.addJob(nowCommand)
			nowCommand = ''
#print len(Jobs.scripts),Jobs.scripts 
Jobs.submit()
