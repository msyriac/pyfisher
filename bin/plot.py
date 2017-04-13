import numpy as np
import sys
import matplotlib.pyplot as plt
import ConfigParser
from makeDerivs import getHubbleCosmology
from glob import glob
#def plotStepSize(paramName):
    #np.loadtxt("output/StepSize/"+str(paramName)+"_"+str(paramVal)+"_"+out_pre+spec+"_lmax_"+str(int(lmax))+".csv",delimiter=",")
def getHubblevsTheta(fparams,step):
    ftheta = fparams['theta100']
    start = ftheta - step
    end   = ftheta + step
    width = step/100.
    ParamRange = np.arange(start,end+width,width)
    data=[]
    for theta in ParamRange:
	H0 = getHubbleCosmology(theta,fparams)
	data.append([theta,H0])
    return np.array(data)


def plotFromFile(filename,xcol,ycol):
    data = np.loadtxt(filename,delimiter=",")
    plt.plot(data[:,xcol],data[:,ycol])
    plt.savefig(filename[0:filename.find('.')]+".png",format='png')
    plt.close()

def getClvsParam(paramName,out_pre,spec,ell):
    fnames=glob("output/StepDump/"+paramName+"_*_"+out_pre+spec+"_lmax_8000.csv")
    data=[]
    for fname in fnames:
	f = np.loadtxt(fname,delimiter=",")
	paramVal = float(fname.split('_')[1])
	data.append(np.insert(f[ell,:],0,paramVal))
    data=np.sort(np.array(data),axis=0)
    fig, ax = plt.subplots(2,3,figsize=(20,10))
    figName = "C_l vs. "+paramName+" (ell="+str(ell)+")"
    fig.suptitle(figName, fontsize=20)
    #TT
    ax[0,0].plot(data[:,0],data[:,1],marker='x')
    ax[0,0].set_title("C_TT")
    #EE
    ax[0,1].plot(data[:,0],data[:,2],marker='x')
    ax[0,1].set_title("C_EE")
    #BB
    ax[0,2].plot(data[:,0],data[:,3],marker='x')
    ax[0,2].set_title("C_BB")
    #TE
    ax[1,0].plot(data[:,0],data[:,4],marker='x')
    ax[1,0].set_title("C_TE")
    #kk
    ax[1,1].plot(data[:,0],data[:,5],marker='x')
    ax[1,1].set_title("C_kk")
    #kt
    ax[1,2].plot(data[:,0],data[:,6],marker='x')
    ax[1,2].set_title("C_kt")
    plt.savefig(figName+out_pre+spec+".png",format='png')
    #plt.show()
    plt.close()
    print "Saved figure ", figName


def main(argv):
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

    #Get data Hubble vs. Theta
    '''
    params = fparams.copy()
    data=getHubblevsTheta(params,stepSizes['theta100'])
    np.savetxt("output/H0vsTheta100/HubbleTheta"+spec+"_lmax"+str(int(fparams['lmax']))+".csv",data,delimiter=",")
    plotFromFile("output/H0vsTheta100/HubbleTheta"+spec+"_lmax"+str(int(fparams['lmax']))+".csv",0,1)
    '''
    #Plot Step Size for theta - input parameter name and ell for plotting!

    paramName='theta100'
    ell = 1000
    for paramName in ['r']:
        for ell in range(0,int(fparams['lmax'])+1,100):
            data = getClvsParam(paramName,out_pre,spec,ell)


if (__name__ == "__main__"):
    main(sys.argv[1:])


