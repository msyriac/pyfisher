import numpy as np
import os,sys
import matplotlib.pyplot as plt
#import matplotlib.lines as mlines
import itertools

testParam = 'w'
ells = np.linspace(0,7000,71)
#ells = [3000]
#ells = np.linspace(2100,3000,10)
#dataRoot = '/gpfs/scratch/nhnguyen/oldtestStepSize/'
dataRoot = os.environ['FISHER_DIR']+'/output/StepDump/'
savename = 'June6_6kmedium_oldtestStepSize_'
name = savename+'highAcc_unlensed_scalar_'

testList = ['w']

#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
#lines = itertools.cycle(['o', 's', '^', 'h', '*', 'd', 'p'])

#xdata = np.linspace(0.0597,0.0603,1001)
#xdata = np.linspace(0.036,0.084,161)
#xdata = np.linspace(-1.4,-0.6,161)
#pmin = -1.4
#xdata = np.linspace(-1.005,-0.995,101)
#pmin = -1.005
xdata = np.linspace(-1.05,-0.95,201)
pmin = -1.05

# round the value to about 6 significant digits                                                                                                          
decimal = int(-(np.ceil(np.log10(abs(pmin)))-6))
xdata = [round(x,decimal) for x in xdata]
print xdata
#sys.exit()

data = np.zeros([len(xdata),7])
#print data
for ell in ells:
    ell = int(ell)
    fig, ax = plt.subplots(2,3,figsize=(20,10))
    figName = "C_l vs. "+testParam+" (ell="+str(ell)+")"
    fig.suptitle(figName, fontsize=20)

    for i in range(len(xdata)):
        ydata = np.loadtxt(dataRoot+name+testParam+'_'+str(xdata[i])+'_fCls.csv',delimiter=',')[ell,:]
        data[i,0] = xdata[i]
        for j in range(6):
            data[i,j+1] = ydata[j]
    #print data
    #sys.exit()
    fileName = os.environ['FISHER_DIR']+"output/StepSizePlots/"+savename+testParam+"_highAcc_ell="+str(ell)
    #np.savetxt(fileName+".csv",data,delimiter=',')
    ax[0,0].plot(data[:,0],data[:,1])
    ax[0,0].set_title("C_TT")
    ax[0,1].plot(data[:,0],data[:,2])
    ax[0,1].set_title("C_EE")
    ax[0,2].plot(data[:,0],data[:,3])
    ax[0,2].set_title("C_BB")
    ax[1,0].plot(data[:,0],data[:,4])
    ax[1,0].set_title("C_TE")
    ax[1,1].plot(data[:,0],data[:,5])
    ax[1,1].set_title("C_kk")
    ax[1,2].plot(data[:,0],data[:,6])
    ax[1,2].set_title("C_kt")

    
    plt.legend(loc='lower left',ncol=3)
    #plt.show()
    plt.savefig(fileName+'.png',format='png')
    print "Saved figure ",fileName
    plt.close()
