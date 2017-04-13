import numpy as np
import sys,os
import matplotlib.pyplot as plt
from makeDerivs import getPowerCamb
import ConfigParser
import itertools

colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y'])

'''
iniFile = os.environ['FISHER_DIR']+'/input/cambRed_test.ini'
Config = ConfigParser.SafeConfigParser()
Config.optionxform = str
Config.read(iniFile)
spec = Config.get('general','spec')
out_pre = Config.get('general','output_prefix')
AccuracyBoost = Config.getboolean('general','AccuracyBoost')
fparams = {}
for (key, val) in Config.items('camb'):
    fparams[key] = float(val)
Cls = getPowerCamb(fparams,spec,AccuracyBoost=AccuracyBoost)
np.savetxt(os.environ['FISHER_DIR']+'/output/test_lensing_pycambCls.dat',Cls)

log = True
prefix = ''
name = 'test1_lensing_scalCovCls'
data = np.loadtxt('output/'+name+'.dat')
ell=1

n = float(len(data[0,:])-ell)
a = int(round(np.sqrt(n),0))
b = int(np.ceil(n/a))
fig, ax = plt.subplots(a,b,figsize=(20,10))
x = data[:,0]
#x = np.arange(float(len(data)))
i = ell

for xpos in range(a):
    for ypos in range(b):
        ax[xpos,ypos].plot(x,data[:,i])
        #ax[xpos,ypos].plot(x,data[:,i]*2*np.pi/(x*(x+1.0)))
        
        if i <= 3:
            ax[xpos,ypos].plot(x,data[:,i]/2/np.pi*(x*(x+1.0)))
        else:
            #ax[xpos,ypos].plot(x,data[:,i]*(x*(x+1.0))**2)
            ax[xpos,ypos].plot(x,data[:,i]*4/2/np.pi)
        
        if log:
            ax[xpos,ypos].set_xscale('log')
            #ax[xpos,ypos].set_yscale('log')
            prefix = '(log)'
        i+=1
        if i>n:
            break
'''
'''
name1 = 'test_lensing_pycambCls'
name2 = 'z1.0_scalCovCls'
name22 = 'z1.0_lenspotentialCls'
name3 = 'z1.0_lensedCls'

data1 = np.loadtxt('output/'+name1+'.dat')
data2 = np.loadtxt('output/'+name2+'.dat')
data22 = np.loadtxt('output/'+name22+'.dat')
data3 = np.loadtxt('output/'+name3+'.dat')

fig = plt.figure()
ax = fig.add_subplot(111)

x1 = np.arange(float(len(data1)))
y1 = data1[:,3]
#y1 = data1[:,4]/x1**4
#y1 = data1[:,5]*np.sqrt(7.4311e12)*2/np.sqrt(x1*(x1+1.))


x2 = data2[:,0]
y2 = data2[:,2]*2.*np.pi/(x2*(x2+1.))
#y2 = data2[:,11]*2.*np.pi/4.*(x2*(x2+1.))
#y2 = data2[:,3]*2.*np.pi/4.*np.sqrt(x2*(x2+1.))
#y2 = data2[:,3]*2.*np.pi/2.

x22 = data22[:,0]
#y22 = data22[:,5]*2*np.pi/(x2*(x2+1))**2
#y22 = data22[:,5]*2.*np.pi/4.
y22 = data22[:,6]*2.*np.pi/4.
y22 = data22[:,6]*2.*np.pi/2./np.sqrt(x22*(x22+1.))

x3 = data3[:,0]
#y3 = data3[:,4]/x3**4/7.4311e12/4
#y3 = data3[:,4]/x3**4/7.4311e12/4*(x3*(x3+1))**2
y3 = data3[:,4]


for label in ['cc','cg','cs','gg','gs','ss']:
    data = np.loadtxt('output/Das_'+label+'.csv',delimiter=',')
    ax.plot(data[:,0],data[:,1],colors.next()+'--',label=label)

names = ['July12_highAcc_Das_bias1_galaxy_fCls']#,'July12_highAcc_Das_bias2_galaxy_fCls']
names = ['July12_highAcc_Das_galaxy_fCls','July12_highAcc_Das_bias1_galaxy_fCls']

lss = itertools.cycle(['-','-.'])
for name in names:
    ls = lss.next()
    data = np.loadtxt('output/'+name+'.csv',delimiter=',')
    x = np.arange(float(len(data)))
    #y = data[:,2]
    
    ax.plot(x,data[:,0],colors.next()+ls,label='$\kappa_{CMB} \kappa_{CMB}$')
    ax.plot(x,data[:,1],colors.next()+ls,label='$\kappa_{CMB} \kappa_{gal}$')
    ax.plot(x,data[:,2],colors.next()+ls,label='$\kappa_{CMB} \Sigma$')
    ax.plot(x,data[:,3],colors.next()+ls,label='$\kappa_{gal} \kappa_{gal}$')
    ax.plot(x,data[:,4],colors.next()+ls,label='$\kappa_{gal} \Sigma$')
    ax.plot(x,data[:,5],colors.next()+ls,label='$\Sigma \Sigma$')
    #ax.plot(x1,y1,label='pycamb')
    #ax.plot(x2,y2,label='scalCov')
    #ax.plot(x22,y22,label='lenspotentialCov')
    #ax.plot(x3,y3,label='lensed')
    #,x2,y2,x3,y3)

ax.set_xscale('log')
ax.set_yscale('log')

#plt.legend(loc='lower left',ncol=2)
plt.show()
fig.suptitle(name+prefix,fontsize=20)
fileName = 'output/'+name+prefix+'.png'
#plt.savefig(fileName,format='png')
#print "Saved figure ",fileName
'''
'''
Cls = np.loadtxt('output/Nov3_highAcc_CDM_unlensed_axion_fCls.csv',delimiter=',')
l = np.array(range(len(Cls)))
imax = 3

for i in range(imax):
    if i == 5:
        Cls[:,i] = Cls[:,i] * np.sqrt(l*(l+1.))
    elif i == 4:
        continue
    else:
        Cls[:,i] = Cls[:,i] * (l*(l+1.))
print Cls

lmin = 0
fig = plt.figure()
for i in range(imax):
    ax = fig.add_subplot(2,3,i+1)
    ax.plot(l[lmin:],Cls[lmin:,i])
    ax.set_xscale('log')
    if i in [4,5]:
        ax.set_yscale('log')
'''
#Cls = np.loadtxt(os.environ['CMBPROJ_DIR']+'/data/TheorySpectra/ell28k_highacc_lensedCls.dat')
Cls = np.loadtxt(os.environ['CMBPROJ_DIR']+'/data/TheorySpectra/Nov10_highAcc_CDM_lensedCls.dat')
#Cls = np.loadtxt(os.environ['AXIONCAMB_DIR']+'Nov8_highAcc_CDM_scalCls.dat')
print Cls
lmin = 6000
fig = plt.figure()
j = 1
for i in range(1,5):
    ax = fig.add_subplot(2,4,j)
    ax.plot(Cls[lmin:,0],Cls[lmin:,i])
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    j+=1
    #vmax = max(Cls[lmin:,i])
    #ax.set_ylim([-vmax/1.e3,vmax])
#Cls = np.loadtxt('/home/nhnguyen/CAMB-May2016/Nov18_highAcc_CDM_lensedCls.dat')    
Cls = np.loadtxt(os.environ['CMBPROJ_DIR']+'/data/TheorySpectra/ell28k_highacc_lensedCls.dat')
#Cls = np.loadtxt(os.environ['CMBPROJ_DIR']+'/data/TheorySpectra/Nov9_highAcc_CDM_lensedCls.dat')
for i in range(1,5):
    ax = fig.add_subplot(2,4,j)
    ax.plot(Cls[lmin:,0],Cls[lmin:,i])
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    j+=1
print Cls
'''
i+=1
Cls = np.loadtxt('output/Nov3_highAcc_CDM_unlensed_axion_fCls.csv',delimiter=',')
x = np.array(range(len(Cls)))
N = len(Cls)
print x
data = x.reshape([N,1])
for j in range(4):
    ax = fig.add_subplot(2,4,i)
    y = Cls[:,j]*(x*(x+1.))/2./np.pi
    ax.plot(x[lmin:],y[lmin:])
    ax.set_xscale('log')
    i+=1
    data = np.hstack([data,y.reshape([N,1])])
    #print y
print data
'''
plt.show()
