import os, sys
import glob
import time
import subprocess
import ConfigParser
import numpy as np
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

paramList = []
stepSizes = {}
fparams = {}
for (key, val) in config.items('camb'):
	if ',' in val:
                param,step = val.split(',')
		paramList.append(key)
		fparams[key]=float(param)
		stepSizes[key]=float(step)

testRange = {}
testList = []

for (key, val) in config.items('testmnu'):
        pmin,pmax,numstep = [float(x) for x in val.split(',')]
	testList.append(key)
	pcenter = fparams[key]
	mstep = (pcenter-pmin)/(numstep/2)
        pstep = (pmax-pcenter)/(numstep/2)
        testRange[key] = np.append(np.arange(pmin,pcenter,mstep),np.arange(pcenter,pmax+pstep/2,pstep))
	# round the value to about 5 significant digits
	decimal = int(-(np.ceil(np.log10(pmin))-5))
        testRange[key] = [round(x,decimal) for x in testRange[key]]
	#print key,pmin,pcenter,pmax

#testRange = {'H0':[67.31]}
#testList = ['H0']	

# a magical formula
numCores=4
N = 0
for testParam in testList:
        # for each value in testRange, needs 1 call for fid and 1 call for deriv for each param
        N+=len(testRange[testParam])*(1+len(paramList))
	print testParam,testRange[testParam]
print "Number of tasks: ", N
tasksPerCore = np.ceil(float(N)/float(numCoresAsked)*8/numCores)
numCoresNeed = float(N) / tasksPerCore
print "Tasks per core: ",tasksPerCore
print "Needed cores: ",int(np.ceil(numCoresNeed)),"/",numCoresAsked
print "Needed nodes: ",int(np.ceil(numCoresNeed/numCores)),"/",numCoresAsked/8
#sys.exit()

commandPreFix = 'source ~/.bashrc'
Jobs = jobMaker(projectName=projectName,commandPreFix=commandPreFix,numCores=numCores,queue=queue,walltime=walltime,jobRoot='/gpfs/scratch/nhnguyen/testmnuDump/jobs/')
optionsCommon = ' --spec ' + spec + ' --out_pre ' + out_pre + ' --AccuracyBoost ' + str(AccuracyBoost)
i = 0
# jobNum = 0
nowCommand = ''
prefix = ''
for testParam in testList:
        for testVal in testRange[testParam]:
                prefix = testParam+'_'+str(testVal)

                # Save fiducials
                i+=1
                saveSuffix = out_pre+spec+"_"+prefix+"_fCls.csv"
                options = optionsCommon + ' --testParam ' + testParam + ' --testVal ' + str(testVal) + ' --paramName 0 --stepSize 0 --output ' + saveRoot+saveSuffix 
                nowCommand += ' python '+os.environ['FISHER_DIR']+'bin/cluster-submit/testmnu_script.py '+options+';'
		if i%tasksPerCore == 0 or i==N:
                        #commandSuff = '"&'
                        #sys.exit()
			Jobs.addJob(nowCommand)
                        nowCommand = ''
                        # jobNum+=1
                # Calculate Derivatives
                for paramName in paramList:
                        i+=1
                        saveSuffix = out_pre+spec+"_"+prefix+"_dCls_"+paramName+".csv"
                        options = optionsCommon+' --testParam '+testParam+' --testVal '+str(testVal)+' --paramName '+paramName+' --stepSize '+str(stepSizes[paramName])+' --output '+saveRoot+saveSuffix 
                        nowCommand += ' python '+os.environ['FISHER_DIR']+'bin/cluster-submit/testmnu_script.py '+options+';'
			if i%tasksPerCore == 0 or i==N:
				Jobs.addJob(nowCommand)
				nowCommand = ''
                                # jobNum+=1
#print len(Jobs.scripts),Jobs.scripts
Jobs.submit()
