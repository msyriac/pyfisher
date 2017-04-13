import numpy as np
import os
import matplotlib.pyplot as plt
import itertools

index = 0
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
#lines = itertools.cycle(['o', 's', '^', 'h', '*', 'd', 'p'])

root = os.environ['FISHER_DIR']+'/output/'
deriv1 = 'June4_optimal_veryhighhighAcc_unlensed_scalar_dCls_mnu.csv'
deriv2 = 'June4_5_testStepSize_veryhighhighAcc_unlensed_scalar_mnu_0.0204_dCls_mnu.csv'
deriv3 = 'June4_3_testStepSize_veryhighhighAcc_unlensed_scalar_mnu_0.000121_dCls_mnu.csv'
deriv5 = 'June4_2_testStepSize_veryhighhighAcc_unlensed_scalar_mnu_2.04e-05_dCls_mnu.csv'
deriv6 = 'June4_3_testStepSize_veryhighhighAcc_unlensed_scalar_mnu_0.000204_dCls_mnu.csv'
deriv7 = 'June4_4_testStepSize_veryhighhighAcc_unlensed_scalar_mnu_0.00204_dCls_mnu.csv'

fid = 'June4_optimal_veryhighhighAcc_unlensed_scalar_fCls.csv'

deriv0 = 'May31_highAcc_unlensed_scalar_dCls_mnu.csv'
fid0 = 'May31_highAcc_unlensed_scalar_fCls.csv'

deriv4 = 'June6_nmn3_vhhAcc_unlensed_scalar_dCls_mnu.csv'
fid4 = 'June6_nmn3_vhhAcc_unlensed_scalar_fCls.csv'


alexdata = np.loadtxt('output/AlexdCls.csv',delimiter=',')
dCls1 = np.loadtxt(root+deriv1,delimiter=',')
dCls2 = np.loadtxt(root+deriv2,delimiter=',')
dCls3 = np.loadtxt(root+deriv3,delimiter=',')
dCls5 = np.loadtxt(root+deriv5,delimiter=',')
dCls6 = np.loadtxt(root+deriv6,delimiter=',')
dCls7 = np.loadtxt(root+deriv7,delimiter=',')
fCls = np.loadtxt(root+fid,delimiter=',')

dCls0 = np.loadtxt(root+deriv0,delimiter=',')
fCls0 = np.loadtxt(root+fid0,delimiter=',')

dCls4 = np.loadtxt(root+deriv4,delimiter=',')
fCls4 = np.loadtxt(root+fid4,delimiter=',')

print dCls1[:3000,index]
print dCls2[:3000,index]
print dCls0[:3000,index]

#print fCls
#print fCls[:3000,index]
#print fCls0[:3000,index]


'''
row = len(fCls)
row = 3000
data = np.zeros([row,6])
for ell in range(row):
    data[ell,0]=ell
    #data[ell,1]=4.0*abs(dCls1[ell,4])/(ell*(ell+1))**2/fCls[ell,4]
    #data[ell,1]=dCls2[ell,1]*2.0*np.pi/(ell*(ell+1))/fCls[ell,1]
    data[ell,1]=dCls1[ell,4]/fCls[ell,4]
    data[ell,2]=dCls2[ell,4]/fCls[ell,4]
    data[ell,3]=dCls3[ell,4]/fCls[ell,4]
    data[ell,4]=dCls0[ell,4]/fCls0[ell,4]

    #data[ell,1] = abs(dCls[ell,4])*2.0/np.pi/fCls[ell,4]
data = np.nan_to_num(data)
#print data
'''
ells,dcls = np.loadtxt('input/cmbpol.csv',unpack=True,delimiter=',')

fig = plt.figure(figsize=(16,6))
ax = fig.add_subplot(111)
#ax.set_title("d", fontsize=16)
#ax.plot(alexdata[:,0],alexdata[:,1],ls='-',label='0.02-Alex')
#ax.plot(ells,dcls/31.0,ls='-',label='cmbpol')

#ax.plot(range(len(fCls)),dCls1[:,index]/fCls[:,index],ls='-',label='4.2e-5')
ax.plot(range(len(fCls)),dCls2[:,index]/fCls[:,index],ls='-',label='0.0204')
#ax.plot(range(len(fCls)),dCls3[:,index]/fCls[:,index],ls='-',label='0.000121')
ax.plot(range(len(fCls)),dCls5[:,index]/fCls[:,index],ls='-',label='2.04e-5')
#ax.plot(range(len(fCls)),dCls6[:,index]/fCls[:,index],ls='-',label='0.000204')
#ax.plot(range(len(fCls)),dCls7[:,index]/fCls[:,index],ls='-',label='0.00204')

#ax.plot(range(len(fCls0)),dCls0[:,index]/fCls0[:,index],ls='-',label='0.02-oldAcc')
#ax.plot(range(len(fCls4)),dCls4[:,index]/fCls4[:,index],ls='-',label='nmn=3')

ax.set_ylabel('$dlnC_{l}^{\phi\phi}/dX$',fontsize=16)
ax.set_xlim([0,3000])
#ax.set_ylim([-0.002,0.002])
#ax.set_xscale('log')
#ax.set_yscale('log')
fileName = root+'June4_dlnClphi.png'
plt.legend(loc='lower right')
plt.show()
#plt.savefig(fileName,format='png')
#print "Saved figure ",fileName

