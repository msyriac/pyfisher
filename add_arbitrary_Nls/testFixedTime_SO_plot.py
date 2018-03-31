import numpy as np
import matplotlib.pyplot as plt
import itertools
color=itertools.cycle(['b','r','g','c','m','y','k'])

#taus = [0.06,0.01,0.005,0.002]
taus = [0.01,0.005,0.002]

prefix = 'output/Feb7_noatm'
prefix_file = 'output/Feb7_file_noatm'


c = color.next()
mnus = np.loadtxt('output/Feb10_atm_tau0.01_mnu.csv')
plt.plot(mnus[:,0],mnus[:,1],c+'-',label='$\\sigma(\\tau)=$0.01 atm')
mnus = np.loadtxt('output/Feb10_atm_file_tau0.01_mnu.csv')
plt.plot(mnus[:,0],mnus[:,1],c+'--',label='$\\sigma(\\tau)=$0.01 atm file')

c = color.next()
mnus = np.loadtxt('../namTest/2uK_atm_mnu.csv')
plt.plot(mnus[:,0],mnus[:,1],c+'o',label='check atm')

for tau in taus:
    c = color.next()
    mnus = np.loadtxt(prefix+'_tau'+str(tau)+'_mnu.csv')
    plt.plot(mnus[:,0],mnus[:,1],c+'-',label='$\\sigma(\\tau)=$'+str(tau))
    mnus = np.loadtxt(prefix_file+'_tau'+str(tau)+'_mnu.csv')
    plt.plot(mnus[:,0],mnus[:,1],c+'--',label='$\\sigma(\\tau)=$'+str(tau)+' file')
c = color.next()
mnus = np.loadtxt('../namTest/2uK_noatm_mnu.csv')
plt.plot(mnus[:,0],mnus[:,1],c+'o',label='check')


plt.legend()
plt.savefig('tests/Feb14_fixedTime_mnu.png')    

plt.clf()

c = color.next()
sns = np.loadtxt('output/Feb10_atm_tau0.01_sn.csv')
plt.plot(sns[:,0],sns[:,1],c+'-',label='$\\sigma(\\tau)=$0.01 atm')
sns = np.loadtxt('output/Feb10_atm_file_tau0.01_sn.csv')
plt.plot(sns[:,0],sns[:,1],c+'--',label='$\\sigma(\\tau)=$0.01 atm file')

c = color.next()
sns = np.loadtxt('../namTest/2uK_atm_sn.csv')
plt.plot(sns[:,0],sns[:,1],c+'o',label='check atm')

for tau in taus:
    c = color.next()
    sns = np.loadtxt(prefix+'_tau'+str(tau)+'_sn.csv')
    plt.plot(sns[:,0],sns[:,1],c+'-',label='$\\sigma(\\tau)=$'+str(tau))
    sns = np.loadtxt(prefix_file+'_tau'+str(tau)+'_sn.csv')
    plt.plot(sns[:,0],sns[:,1],c+'--',label='$\\sigma(\\tau)=$'+str(tau)+' file')
c = color.next()
sns = np.loadtxt('../namTest/2uK_noatm_sn.csv')
plt.plot(sns[:,0],sns[:,1],c+'o',label='check')



plt.legend()
plt.savefig('tests/Feb14_fixedTime_sn.png')    

#plt.show()
#sns = np.loadtxt('Jan29_tau'+str(tau)+'_noatm_sn.csv')
