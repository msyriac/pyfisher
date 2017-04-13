import os
import itertools
import numpy as np
import sys
import ConfigParser
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import matplotlib

def confEllipse(param1,param2,sigmax2,sigmay2,sigmaxy,xcenter,ycenter,confLevel=1,savefig=True,savedata=False,verbose=False):
    script = '[confEllipse]'
    alpha = {1:1.52, 2:2.48, 3:3.41}
    # Ellipse parameters                                                                                                                                     
    if ((sigmay2/sigmax2)<1e-10) or ((sigmax2/sigmay2)<1e-10):
        a = alpha[confLevel]*np.sqrt(max(sigmax2,sigmay2) + sigmaxy**2/max(sigmax2,sigmay2))
        b = alpha[confLevel]*np.sqrt(min(sigmax2,sigmay2) - sigmaxy**2/max(sigmax2,sigmay2))
    else:
        a = alpha[confLevel]*np.sqrt((sigmax2+sigmay2)/2. + np.sqrt( (sigmax2-sigmay2)**2/4. + sigmaxy**2 ))
        b = alpha[confLevel]*np.sqrt((sigmax2+sigmay2)/2. - np.sqrt( (sigmax2-sigmay2)**2/4. + sigmaxy**2 ))
    angle = 1./2.*np.arctan(2.*sigmaxy/(sigmax2-sigmay2))*180./np.pi
    if (sigmax2<sigmay2):
        c = a
        a = b
        b = c
    width = 2.*a
    height = 2.*b
    
    script += '\nsigmax2 = '+str(sigmax2)
    script += '\nsigmay2 = '+str(sigmay2)
    script += '\nsigmaxy = '+str(sigmaxy)
    script += '\nparams = '+param1+','+param2
    script += '\nconfLevel = '+str(confLevel)
    script += '\nangle = '+str(angle)
    script += '\nxcenter = '+str(xcenter)
    script += '\nycenter = '+str(ycenter)
    script += '\nwidth = '+str(width)
    script += '\nheight = '+str(height)
    if verbose:
        print(script)
    fileName = os.environ['FISHER_DIR']+'/output/July6_confEllipse_'+param1+'_'+param2+'_'+str(confLevel)+'sigma_delens'
    
    if savedata:
        with open(fileName+'.csv','w') as tempFile:
            tempFile.write(script)
        print('Saved file '+fileName+'.csv')
'''
param1 = 'ns'
param2 = 'r'
sigmax2 = 0.002**2
sigmay2 = (0.95e-3)**2
sigmaxy = 0.0
xcenter = 0.9655
ycenter = 0.0
confEllipse(param1,param2,sigmax2,sigmay2,sigmaxy,xcenter,ycenter,confLevel=1,savedata=True,verbose=True)
sys.exit()
'''
colors = itertools.cycle(['b', 'r', 'g', 'm', 'y', 'c', 'k'])
matplotlib.rcParams['mathtext.default'] = 'regular'
labels = {'H0':'$H_0$','ombh2':'$\Omega_b h^2$','omch2':'$\Omega_c h^2$','ns':'$n_s$','As':'$A_s$','tau':'$\\tau$','mnu':'$\Sigma m_{\\nu}$','nnu':'$N_{eff}$','r':'$r$'}
CL = {1:'68%',2:'95%',3:'99%'}
fontsize = 16

#dataFile1 = os.environ['FISHER_DIR']+'/output/July6_confEllipse_ns_r_1sigma_nodelens.csv'
#dataFile2 = os.environ['FISHER_DIR']+'/output/July6_confEllipse_ns_r_1sigma_delens.csv'
dataFile1 = os.environ['FISHER_DIR']+'/output/June7_newoptimal_vhhAcc_unlensed_scalar_Planck4pt_confEllipse_mnu_nnu_1sigma.csv'
dataFile2 = os.environ['FISHER_DIR']+'/output/June7_newoptimal_vhhAcc_lensed_scalar_S4Primary_confEllipse_mnu_nnu_1sigma.csv'
dataFile3 = os.environ['FISHER_DIR']+'/output/June7_newoptimal_vhhAcc_unlensed_scalar_S4Primary+S4Lens_confEllipse_mnu_nnu_1sigma.csv'
dataFile4 = os.environ['FISHER_DIR']+'/output/June7_newoptimal_vhhAcc_lensed_scalar_S4Primary+DESI_confEllipse_mnu_nnu_1sigma.csv'
dataFile5 = os.environ['FISHER_DIR']+'/output/June7_newoptimal_vhhAcc_unlensed_scalar_S4Primary+S4Lens+DESI_confEllipse_mnu_nnu_1sigma.csv'

#dataFile6 = os.environ['FISHER_DIR']+'/output/Sep24_vhhAcc_unfixKT_unlensed_axion_S4Primary+S4Lens_confEllipse_mnu_nnu_1sigma.csv'
dataFile6 = os.environ['FISHER_DIR']+'/output/Sep28_vhhAcc_pyCAMBunlensed_scalar_CMB_confEllipse_mnu_nnu_1sigma.csv'
dataFile7 = os.environ['FISHER_DIR']+'/output/Sep28_vhhAcc_pyCAMBunlensed(fix)_scalar_CMB_confEllipse_mnu_nnu_1sigma.csv'
#dataFile7 = os.environ['FISHER_DIR']+'/output/Sep28_vhhAcc_unfixKT_unlensed_axion_CMB_confEllipse_mnu_nnu_1sigma.csv'
#dataFile8 = os.environ['FISHER_DIR']+'/output/Sep28_vhhAcc_fixKT_unlensed_axion_CMB_confEllipse_mnu_nnu_1sigma.csv'


dataFiles = [dataFile1,dataFile2,dataFile3,dataFile4,dataFile5,dataFile6,dataFile7]

dataLabels = ['Planck 2pt + Planck Lens 4pt','S4 2pt','S4 2pt + S4 Lens 4pt','S4 2pt + DESI','S4 2pt + S4 Lens 4pt + DESI','py unfix','py fix'] #,'axion unfix sep24','axion unfix sep28','axion fix sep28']
#dataLabels = ['S4 with no delensing','S4 with delensing']
alpha = {1:1.52, 2:2.48, 3:3.41}

#fig = plt.figure(figsize=(10,6))
fig = plt.figure()
ax = fig.add_subplot(111)
'''
#Das3 = np.loadtxt('output/Das_CMB_omch2_w_param.csv',delimiter=',')
#Das4 = np.loadtxt('output/Das_CMB+Lensr_omch2_w_param.csv',delimiter=',')
Das5 = np.loadtxt('output/Das_omL_w_param.csv',delimiter=',')
Das6 = np.loadtxt('output/Das_LensR_omL_w_param.csv',delimiter=',')
#Das7 = np.loadtxt('output/Das_LensE_omL_w_param.csv',delimiter=',')
#Das8 = np.loadtxt('output/Das_LensR+LensE_omL_w_param.csv',delimiter=',')

DasLabels=['Planck CMB','Planck CMB+1% Lens ratio']
Das = [Das5,Das6]
#dxdy=0.01959886929600215
dydx=2.7322695035460995
for i in range(len(Das)):
    data = Das[i]
    dataLabel = DasLabels[i]
    height = np.sqrt(((data[0,0]-data[1,0]))**2+((data[0,1]-data[1,1]))**2)
    width  = np.sqrt(((data[2,0]-data[3,0]))**2+((data[2,1]-data[3,1])/dydx)**2)
    angle = np.arctan((data[2,1]-data[3,1])/dydx**2/(data[2,0]-data[3,0]))*180.0/np.pi
    #angle = np.arctan((data[1,1]-data[0,1])/dydx/(data[1,0]-data[0,0]))*180.0/np.pi
    #angle = np.arctan((data[1,1]-data[0,1])*dxdy/(data[1,0]-data[0,0]))*180.0/np.pi
    print width,height,angle
    xcenter = 0.74
    ycenter = -1.0
   
    for angle in np.linspace(0.,180.,25):
        e = Ellipse((xcenter,ycenter),width,height, angle=angle,color = colors.next(),fill=False,label=dataLabel) 
        ax.add_patch(e)
    
    e = Ellipse((xcenter,ycenter),width,height, angle=angle,color = colors.next(),fill=False,label=dataLabel) 
    ax.add_patch(e)
    
    ax.plot(xcenter,ycenter,'r*')#,markersize=16)
ax.set_xlabel('$\Omega_{\Lambda}$',fontsize=16)
ax.set_ylabel('w',fontsize=16)

'''
#param1 = 'ns'
#param2 = 'r'
 
for i in range(len(dataFiles)):
    dataFile = dataFiles[i]
    dataLabel = dataLabels[i]

    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(dataFile)
  
    sigmax2 = config.getfloat('confEllipse','sigmax2')
    sigmay2 = config.getfloat('confEllipse','sigmay2')
    sigmaxy = config.getfloat('confEllipse','sigmaxy')
    try:
        param1,param2 = config.get('confEllipse','params').split(',')
    except:
        param1 = 'mnu'
        param2 = 'nnu'
    confLevel = config.getint('confEllipse','confLevel')
    angle = config.getfloat('confEllipse','angle')
    xcenter = config.getfloat('confEllipse','xcenter')
    ycenter = config.getfloat('confEllipse','ycenter')
    width = config.getfloat('confEllipse','width')
    height = config.getfloat('confEllipse','height')

    print xcenter,ycenter
    #xcenter = 0.09921434062740242
    #ycenter = 3.0389654547386447
    #xcenter = 0.06
    #ycenter = 0.12
    #ycenter = 0.08
    e = Ellipse((xcenter,ycenter),width,height, angle=angle,color = colors.next(),fill=False,label=dataLabel)
    ax.add_patch(e)
    ax.plot(xcenter,ycenter,'r*')#,markersize=16)

    ax.set_xlabel(labels[param1],fontsize=fontsize)
    ax.set_ylabel(labels[param2],fontsize=fontsize)

#ax.set_yscale('log')
#ax.set_ylim([0,4e-3])
#ax.set_xlim([0.96,0.975])
ax.set_xlim([0,0.4])
ax.set_ylim([2.8,3.3])
plt.grid()
plt.legend(loc='upper right')
ax.set_title('Joint constraint ('+CL[confLevel]+' CL) on '+labels[param1]+' and '+labels[param2],fontsize=fontsize)
fileName = os.environ['FISHER_DIR']+'/output/July6_confEllipse_'+param1+'_'+param2+'_'+str(confLevel)+'sigma'
#fileName = os.environ['FISHER_DIR']+'/output/June29_Das_confEllipse_omL_w_1sigma'
plt.show()
#plt.savefig('Sep28_fixKT',format='png')
#plt.savefig(fileName+'.png',format='png')

#print('Saved file '+fileName+'.png')
            
