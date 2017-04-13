import os, sys
import glob
import time
import subprocess
import ConfigParser
import numpy as np
from gpcInterface import jobMaker
    
iniFile = os.environ['FISHER_DIR']+'/input/cluster-submit/makeDerivs_input.ini'
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

paramList = []
stepSizes = {}
fparams = {}
for (key, val) in config.items('camb'):
	if ',' in val:
                param,step = val.split(',')
		paramList.append(key)
		fparams[key]=float(param)
		stepSizes[key]=float(step)
# a magical formula
numCores=1
N = 0
# for each value in testRange, needs 1 call for fid and 1 call for deriv for each param
N=1+len(paramList)

print "Number of tasks: ", N
tasksPerCore = np.ceil(float(N)/float(numCoresAsked)*8/numCores)
numCoresNeed = float(N) / tasksPerCore
print "Tasks per core: ",tasksPerCore
print "Needed cores: ",int(np.ceil(numCoresNeed)),"/",numCoresAsked
print "Needed nodes: ",int(np.ceil(numCoresNeed/numCores)),"/",numCoresAsked/8
#sys.exit()

commandPreFix = 'source ~/.bashrc'
Jobs = jobMaker(projectName=projectName,commandPreFix=commandPreFix,numCores=numCores,queue=queue,walltime=walltime,jobRoot='/gpfs/scratch/nhnguyen/makeDerivs/jobs/')
optionsCommon = ' --spec ' + spec + ' --out_pre ' + out_pre + ' --AccuracyBoost ' + str(AccuracyBoost)
i = 0
nowCommand = ''

# Save fiducials
i+=1
saveSuffix = out_pre+'_'+spec
options = optionsCommon + ' --paramName 0 --stepSize 0 --output ' + saveRoot+saveSuffix 
nowCommand += ' python '+os.environ['FISHER_DIR']+'bin/cluster-submit/makeDerivs_script.py '+options+';'
if i%tasksPerCore == 0 or i==N:
	Jobs.addJob(nowCommand)
	nowCommand = ''
# Calculate Derivatives
for paramName in paramList:
	i+=1
	#saveSuffix = out_pre+spec+"_dCls_"+paramName+".csv"
	options = optionsCommon+' --paramName '+paramName+' --stepSize '+str(stepSizes[paramName])+' --output '+saveRoot+saveSuffix 
	nowCommand += ' python '+os.environ['FISHER_DIR']+'bin/cluster-submit/makeDerivs_script.py '+options+';'
	if i%tasksPerCore == 0 or i==N:
		Jobs.addJob(nowCommand)
		nowCommand = ''
#print len(Jobs.scripts),Jobs.scripts
Jobs.submit()
