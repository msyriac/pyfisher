import numpy as np
import matplotlib.pyplot as plt

#taus = [0.06,0.01,0.005,0.002]
taus = [0.01,0.005,0.002]

for tau in taus:
    mnus = np.loadtxt('Jan29_tau'+str(tau)+'_noatm_mnu.csv')
    plt.plot(mnus[:,0],mnus[:,1],label='$\\sigma(\\tau)=$'+str(tau))
mnus = np.loadtxt('2uK_noatm_mnu.csv')
plt.plot(mnus[:,0],mnus[:,1],'o',label='check')
mnus = np.loadtxt('Byeonghee_email.csv')
plt.plot(mnus[:,0],mnus[:,2],'o',label='email check')
plt.legend()
plt.savefig('fixedTime_mnu.png')    

plt.clf()
for tau in taus:
    sns = np.loadtxt('Jan29_tau'+str(tau)+'_noatm_sn.csv')
    plt.plot(sns[:,0],sns[:,1],label='$\\sigma(\\tau)=$'+str(tau))
sns = np.loadtxt('2uK_noatm_sn.csv')
plt.plot(sns[:,0],sns[:,1],'o',label='check')
plt.legend()
plt.savefig('fixedTime_sn.png')    

#plt.show()
#sns = np.loadtxt('Jan29_tau'+str(tau)+'_noatm_sn.csv')
