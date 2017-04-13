import numpy as np
import os,sys
import matplotlib.pyplot as plt
#import matplotlib.lines as mlines
import itertools

#testParam = 'mnu'
#testParam = 'H0'
'''
dataRoot1 = os.environ['FISHER_DIR']+'/output/June2_small_testStepSize_highAcc_unlensed_scalar_'
dataRoot4 = os.environ['FISHER_DIR']+'/output/June2_large_testStepSize_highAcc_unlensed_scalar_'
dataRoot3 = os.environ['FISHER_DIR']+'/output/May31_testStepSize_highAcc_unlensed_scalar_'
dataRoot2 = os.environ['FISHER_DIR']+'/output/June1_testStepSize_highAcc_unlensed_scalar_'
'''
#dataRoot = os.environ['FISHER_DIR']+'/output/June3_testStepSize_highAcc_unlensed_scalar_'
#dataRoot = os.environ['FISHER_DIR']+'/output/June4_testStepSize_highhighAcc_unlensed_scalar_'
#dataRoot = os.environ['FISHER_DIR']+'/output/June4_testStepSize_veryhighAcc_unlensed_scalar_'
dataRoot = os.environ['FISHER_DIR']+'/output/June23_testStepSize_240meV_unlensed_scalar_'

#dataRoot1 = os.environ['FISHER_DIR']+'/output/June4_testStepSize_superhighAcc_unlensed_scalar_'
#dataRoot2 = os.environ['FISHER_DIR']+'/output/June4_large_testStepSize_superhighAcc_unlensed_scalar_'

paramList = ['H0','ombh2','omch2','tau','ns','As','mnu']
#paramList = ['H0','ombh2','omch2','tau','ns','As','w']
#testList = ['H0','ombh2','omch2','tau','ns','As','mnu','nnu']
testList = ['mnu']
model = 'LCDM+mnu'
index = 6
paramName = paramList[index]

colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
lines = itertools.cycle(['o', 's', '^', 'h', '*', 'd', 'p','v'])

fig = plt.figure(figsize=(16,6))
#fig = plt.figure(figsize=(16,14))
ax = fig.add_subplot(111)
#ax.set_title("Constraint on neutrino mass from CMB-S4 vs.\nStep size of LCDM+mnu+nnu parameters", fontsize=16)
#ax.set_title("Constraint on equation-of-state parameter from CMB-S4 vs.\nStep size of LCDM+w parameters", fontsize=16)
ax.set_title("Constraint on "+paramName+" from CMB-S4 vs.\nStep size of "+model+" parameters", fontsize=16)
'''
xdata1 = np.logspace(np.log10(0.002),np.log10(0.2),100)
xdata4 = np.arange(31.0,100.5,1.0)
xdata3 = np.arange(1.0,30.5,1.0)
xdata2 = np.arange(0.05,0.97,0.05)
xdata = np.append(np.append(np.append(xdata1,xdata2),xdata3),xdata4)

xdata1 = np.logspace(np.log10(0.001),np.log10(0.01),50)
xdata2 = np.logspace(np.log10(0.01),np.log10(0.1),50)
xdata3 = np.logspace(np.log10(0.1),np.log10(1.0),50)
xdata4 = np.logspace(np.log10(1.0),np.log10(10.0),50)
xdata5 = np.logspace(np.log10(10.0),np.log10(100.0),50)
xdata = np.append(np.append(np.append(np.append(xdata1,xdata2),xdata3),xdata4),xdata5)
'''
xdata = np.logspace(np.log10(0.01),np.log10(100.0),160)
for testParam in testList:
    '''
    ydata1 = np.loadtxt(dataRoot1+testParam+'_'+model+'.csv',delimiter=',')[:,index+1]
    ydata2 = np.loadtxt(dataRoot2+testParam+'_'+model+'.csv',delimiter=',')[:,index+1]
    ydata3 = np.loadtxt(dataRoot3+testParam+'_'+model+'.csv',delimiter=',')[:,index+1]
    ydata4 = np.loadtxt(dataRoot4+testParam+'_'+model+'.csv',delimiter=',')[:,index+1]
    ydata = np.append(np.append(np.append(ydata1,ydata2),ydata3),ydata4)
    '''
    #ydata1 = np.loadtxt(dataRoot1+testParam+'_'+model+'_optimal.csv',delimiter=',')[:,index+1]
    #ydata2 = np.loadtxt(dataRoot2+testParam+'_'+model+'_optimal.csv',delimiter=',')[:,index+1]
    #ydata = np.append(ydata1,ydata2)
    ydata = np.loadtxt(dataRoot+testParam+'_'+model+'_S4.csv',delimiter=',')[:,index+1]
    if paramName == 'mnu':
        ydata = ydata*1e3
    
    data = np.array([xdata,ydata]).T
    data = data[np.lexsort(np.fliplr(data).T)]
    ax.plot(data[:,0],data[:,1],colors.next()+lines.next(),ls='-',label=testParam)
    #ax.plot(xdata,ydata,colors.next()+lines.next(),ls='-',label=testParam)
plt.legend(loc='lower right',ncol=3)
#ax.set_ylim([0.,0.05])
#ax.set_ylim([0,0.210])
#ax.set_ylim([0,0.015])
#ax.set_xlim([1e-3,1e1])
#ax.set_xlim([0.95,30.5])
ax.set_xscale('log')
#ax.set_yscale('log')

#ax.set_ylabel('$\sigma_{H0}$ (km/s/Mpc)',fontsize=16)
#ax.set_ylabel('$\sigma_{mnu}$ (meV)',fontsize=16)
ax.set_ylabel('$\sigma$('+paramName+')',fontsize=16)
#ax.set_xlabel('$\\theta_i$ (LCDM+mnu) step size (% of CMB-S4 best fit)',fontsize=16)
#ax.set_ylabel('$\sigma(\omega)$',fontsize=16)
ax.set_xlabel('$\\theta_i$ ('+model+') step size (% of CMB-S4 best fit)',fontsize=16)
'''
xticks=ax.get_xticks().tolist()
for i in range(len(xticks)):
    if (abs(xticks[i])<1.0e-5):
        xticks[i]='CMB-S4 fid'
    else:
        xticks[i]=str(int(xticks[i]))+'$\sigma$'
ax.set_xticklabels(xticks)
'''
fileName = os.environ['FISHER_DIR']+"output/June23_testStepSize_240meV_"+paramName+"_"+model+"_S4(log).png"
#plt.show()
plt.savefig(fileName,format='png')
print "Saved figure ",fileName
