import sys
import numpy as np
from driver import FisherForecast
import ConfigParser
import matplotlib.pyplot as plt
import os

def main(argv):
    try:
        iniFile = argv[0]
    except:
        iniFile = os.environ['FISHER_DIR']+"/input/testmnu_driver.ini"
	finiFile = os.environ['FISHER_DIR']+'/input/cluster-submit/testmnu_input.ini'

    # Read Config
    fConfig = ConfigParser.SafeConfigParser()
    fConfig.optionxform = str
    fConfig.read(finiFile)


    #paramList = []
    #stepSizes = {}
    fparams = {}
    for (key, val) in fConfig.items('camb'):
	    if ',' in val:
		    param,step = val.split(',')
		    #paramList.append(key)
		    fparams[key]=float(param)
		    #stepSizes[key]=float(step)

    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)

    testRange = {}
    testList = []
    for (key, val) in Config.items('testmnu'):
	    if key == 'include': continue    
	    pmin,pmax,numstep = [float(x) for x in val.split(',')]
	    testList.append(key)
	    pcenter = fparams[key]
	    mstep = (pcenter-pmin)/(numstep/2)
	    pstep = (pmax-pcenter)/(numstep/2)
	    testRange[key] = np.append(np.arange(pmin,pcenter,mstep),np.arange(pcenter,pmax+pstep/2,pstep))
	    # round the value to about 5 significant digits
	    decimal = int(-(np.ceil(np.log10(pmin))-5))
	    testRange[key] = [round(x,decimal) for x in testRange[key]]
	    print key,testRange[key]

    for testParam in testList:
	    print "Testing sigma vs. fiducial ",testParam
	    mnudata = np.zeros([len(testRange[testParam]),2])
	    H0data = np.zeros([len(testRange[testParam]),2])
	    i=0
	    for testVal in testRange[testParam]:
		    print testParam," = ",testVal
		    '''
		    if (abs(fraction-1.0)<1e-10):
		    prefix = '100.0'
		    else:
		    prefix = str(fraction*100)+testParam
		    '''
		    prefix = testParam+'_'+str(testVal)
		    F = FisherForecast(iniFile,prefix)
		    F.calcFisher(verbose = True)
		    sigmamnu = 1.e3*F.margSigma("mnu")
		    sigmaH0 = F.margSigma("H0")
            #print "sigma_mnu = ",'{:3.0f}'.format(sigma), " meV for ",testParam," = ",testVal
		    mnudata[i,0] = testVal
		    H0data[i,0] = testVal
		    mnudata[i,1] = sigmamnu
		    H0data[i,1] = sigmaH0
		    i+=1
	    np.savetxt(os.environ['FISHER_DIR']+"/output/June3_testfid_mnu_"+testParam+".csv",mnudata,delimiter=",")
	    np.savetxt(os.environ['FISHER_DIR']+"/output/June3_testfid_H0_"+testParam+".csv",H0data,delimiter=",")
	    '''
	    # Plot results
	    fig = plt.figure()
	    ax = fig.add_subplot(111)
	    figName = "testmnu_May26_"+testParam
	    ax.set_title(figName, fontsize=20)
	    ax.plot(plotdata[:,0],plotdata[:,1],'*')
	    plt.savefig("output/"+figName+".png",format='png')
	    #plt.show()                                                                                                                        
	    plt.close()
	    print "Saved figure ", figName
	    '''
if (__name__ == "__main__"):
	main(sys.argv[1:])

