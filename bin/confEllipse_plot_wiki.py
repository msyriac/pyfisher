import os,sys
import itertools
import numpy as np
import ConfigParser
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import matplotlib
import argparse

# Pass multiple dataFiles and dataLabels (separated by commas) to plot different confidence ellipses on the same plot
# Parse Arguments
parser = argparse.ArgumentParser(description='Plot Confidence Ellipse',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dataFiles',type=str)
parser.add_argument('--dataLabels',type=str)
args = parser.parse_args()
dataFiles = args.dataFiles.split(',')
dataLabels = args.dataLabels.split(',')

# Hard-coded
matplotlib.rcParams['mathtext.default'] = 'regular'
colors = itertools.cycle(['b', 'g', 'r', 'm', 'y', 'c', 'k'])
labels = {'H0':'$H_0$','ombh2':'$\Omega_b h^2$','omch2':'$\Omega_c h^2$','ns':'$n_s$','As':'$A_s$','tau':'$\\tau$','mnu':'$\Sigma m_{\\nu}$','nnu':'$N_{eff}$','r':'$r$','w':'$w$'}
CL = {1:'68%',2:'95%',3:'99.7%'}
fontsize = 16

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(len(dataFiles)):
    dataFile = dataFiles[i]
    dataLabel = dataLabels[i]

    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(dataFile)
  
    param1,param2 = config.get('confEllipse','params').split(',')
    confLevel = config.getint('confEllipse','confLevel')
    angle = config.getfloat('confEllipse','angle')
    xcenter = config.getfloat('confEllipse','xcenter')
    ycenter = config.getfloat('confEllipse','ycenter')
    width = config.getfloat('confEllipse','width')
    height = config.getfloat('confEllipse','height')

    e = Ellipse((xcenter,ycenter),width,height, angle=angle,color = colors.next(),fill=False,label=dataLabel)
    ax.add_patch(e)
    ax.plot(xcenter,ycenter,'r*')

    ax.set_xlabel(labels[param1],fontsize=fontsize)
    ax.set_ylabel(labels[param2],fontsize=fontsize)

# Hard coded again
ax.set_xlim([0.1185,0.1210])
ax.set_ylim([-1.15,-0.85])

ax.set_title('Joint constraint ('+CL[confLevel]+' CL) on '+labels[param1]+' and '+labels[param2],fontsize=fontsize)
plt.grid()
plt.legend(loc='upper right')
fileName ='output/July14_confEllipse_'+param1+'_'+param2+'_'+str(confLevel)+'sigma'
plt.savefig(fileName+'.png',format='png')
print('Saved file '+fileName+'.png')
#plt.show()
            
