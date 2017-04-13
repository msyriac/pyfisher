import numpy as np
import os,sys
import matplotlib.pyplot as plt
#import matplotlib.lines as mlines
import itertools

#testParam = 'mnu'
#testParam = 'H0'

#dataRoot = os.environ['FISHER_DIR']+'output/June6_testfid_veryhhAcc_unlensed_scalar_'
dataRoot = os.environ['FISHER_DIR']+'output/June6_testfid_vhhAcc_unlensed_scalar_'
paramList = ['H0','ombh2','omch2','tau','ns','As','mnu']
testList = ['H0','ombh2','omch2','tau','ns','As']
#testList = ['mnu']                                                                                                                                              
model = 'LCDM+mnu'
index = 6
paramName = paramList[index]

colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
lines = itertools.cycle(['o', 's', '^', 'h', '*', 'd', 'p'])

fig = plt.figure(figsize=(16,6))
ax = fig.add_subplot(111)
ax.set_title("Constraint on "+paramName+" from CMB-S4 vs.\nFiducial value of LCDM parameters", fontsize=16)
#ax.set_title("Constraint on Hubble constant from CMB-S4 vs.\nFiducial value of LCDM parameters", fontsize=16)

xdata = np.arange(-2.0,2.1,0.2)
for testParam in testList:
    ydata = np.loadtxt(dataRoot+testParam+'.csv',delimiter=',')[:,index+1]
    #print xdata,ydata
    if paramName == 'mnu':
        ydata = ydata*1e3
    ax.plot(xdata,ydata,colors.next()+lines.next(),ls='-',label=testParam)
plt.legend(loc='lower left',ncol=3)
#ax.set_ylim([60.,80.])
#ax.set_ylim([0.7,1.0])
ax.set_xlim([-2.2,2.2])

#ax.set_ylabel('$\sigma_{H0}$ (km/s/Mpc)',fontsize=16)
ax.set_ylabel('$\sigma_{mnu}$ (meV)',fontsize=16)
ax.set_xlabel('$\\theta_i$ (Lambda-CDM)',fontsize=16)

xticks=ax.get_xticks().tolist()
for i in range(len(xticks)):
    if (abs(xticks[i])<1.0e-5):
        xticks[i]='CMB-S4 fid'
    else:
        xticks[i]=str(int(xticks[i]))+'$\sigma$'
ax.set_xticklabels(xticks)
fileName = os.environ['FISHER_DIR']+'output/June6_testfid_'+paramName+'_'+model+'_CMB(0.2%).png'
#plt.show()
plt.savefig(fileName,format='png')
print "Saved figure ",fileName
