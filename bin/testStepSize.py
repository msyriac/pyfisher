import camb
from camb import model, initialpower
import numpy as np
import sys
import ConfigParser
import matplotlib.pyplot as plt
from timeit import default_timer as timer

from makeDerivs import getHubbleCosmology,getPowerCamb

#from mmUtils import Plotter

'''
Needs to make sure fiducial parameters don't change after looping
Remember to check the steps in the papers.
Need to fix str(paramName) to paramName only
'''

#Step Size Test
# ParamRange can be manually created
def StpSzTest(params,out_pre,spec,AccuracyBoost,ParamRange,paramName,lmax,verbose=True):
    #StpSzs=[scale*float(stepSizes[paramName]) for scale in ScaleSet]
    #data=np.zeros((NumStep+1,6))
    #Return list containing "StpSz_x": [parval,TT,EE,BB,TE,kk] that can be plotted
    #result={}
    #Record fiducial value
    #parfidval=params[paramName]
    for paramVal in ParamRange:
        if verbose: print paramName,paramVal
        params[paramName]=paramVal
    	if paramName=='theta100':
	    params['H0'] = getHubbleCosmology(params['theta100'],params)
        CLS = getPowerCamb(params,spec,AccuracyBoost=AccuracyBoost)
	#Name gives value of the parameter and lmax
	np.savetxt("output/StepSize/"+paramName+"_"+str(paramVal)+"_"+out_pre+spec+"_lmax_"+str(int(lmax))+".csv",CLS,delimiter=",")
    print "--------------Done with ",paramName

def main(argv):
    starttime = timer()
    # Read Config
    iniFile = "input/StepSize.ini"
    Config = ConfigParser.SafeConfigParser()
    # Reserve case when get string
    Config.optionxform = str
    Config.read(iniFile)
    spec = Config.get('general','spec')
    out_pre = Config.get('general','output_prefix')
    AccuracyBoost = Config.getboolean('general','AccuracyBoost')
    paramList = []
    fparams = {}
    stepSizes = {}
    testSteps={}
    for (key, val) in Config.items('camb'):
        if ',' in val:
            param, step, testStep = val.split(',')
            paramList.append(key)
            fparams[key] = float(param)
            stepSizes[key] = float(step)
	    testSteps[key] = bool(int(testStep))
        else:
            fparams[key] = float(val) 

    #start,end,width = [float(x) for x in Config.get('StpSzTest','ParamRange').split(',')]

    if not('H0' in fparams):
        fparams['H0'] = getHubbleCosmology(fparams['theta100'],fparams)
    print paramList
    print fparams
    for paramName in paramList:
	if not testSteps[paramName]: continue
	#Make range for each parameter
	#----------------------------
	#First test: x2 the range, x0.01 the stepsize
	start = fparams[paramName] - stepSizes[paramName]
	end   = fparams[paramName] + stepSizes[paramName]
	width = stepSizes[paramName]/100.
	ParamRange = np.arange(start,end+width,width)
	#----------------------------
	StpSzTest(fparams,out_pre,spec,AccuracyBoost,ParamRange,paramName,fparams['lmax'],verbose=True)
	endtime = timer()
	print "Time elapsed: ",endtime-starttime

    # pl = Plotter()
    # pl.add(ParamRange,DataSet[:,0])
    # pl.done("output/params.png")
    
    
    # figSz = int(bool(np.sqrt(len(ScaleSet))%1))
    # fig, ax = plt.subplots(figSz,figSz)
    # fig.suptitle('TT', fontsize=20)
    # ikey=0
    # if verbose: print "Plotting"
    # for i in range(0,figSz-1):
    #     for j in range(0,figSz-1):
    #         data=DataSet[DataSet.keys()[ikey]]
    #         ax[i,j].plot(data[:,0],data[:,1])
    #         ax[i,j].set_title(DataSet.keys()[ikey])
    #         ikey+=1
    #         if verbose: print "Final steps..."
    #         if ikey>len(ScaleSet):
    #             break
    # plt.show() 
    # end = timer()
    # if verbose: print "Time elapsed: ",end-start
if (__name__ == "__main__"):
    main(sys.argv[1:])
