
import numpy as np
import sys,os
import ConfigParser
import itertools
import matplotlib.pyplot as plt

'''
paramNames = ['H0','ombh2','omch2','As','ns','tau','mnu','nnu']
powers = ['TT','EE','BB','TE','KK','KT']

for i in range(len(powers)):
    fig = plt.figure(figsize=(16,10))
    for index in range(len(paramNames)):
        ax = fig.add_subplot(3,3,index+1)
        param = paramNames[index]
        
        axion = np.loadtxt('output/Sep28_vhhAcc_fixKT_unlensed_axion_dCls_'+param+'.csv',delimiter=',')
        #axion = np.loadtxt('output/Sep28_vhhAcc_pyCAMBunlensed(fix)_scalar_dCls_'+param+'.csv',delimiter=',')
        python = np.loadtxt('output/Sep28_vhhAcc_unfixKT_unlensed_axion_dCls_'+param+'.csv',delimiter=',')
        #python = np.loadtxt('output/Sep28_vhhAcc_pyCAMBunlensed_scalar_dCls_'+param+'.csv',delimiter=',')
        #python= np.loadtxt('output/June7_newoptimal_vhhAcc_unlensed_scalar_dCls_'+param+'.csv',delimiter=',')
        i = 5
        x1 = np.array(range(len(axion[:,i])))
        y1 = axion[:,i]
        x2 = np.array(range(len(python[:,i])))
        y2 = python[:,i]
        if i == 5:
            #y2 = np.nan_to_num(y2 / np.sqrt(x2*(x2+1.)) *(2.7255*1.e6) * 2)
            y1 = np.nan_to_num(y1 * np.sqrt(x1*(x1+1.)) /(2.7255*1.e6) / 2)
            #python[:,i] = y2
            #axion[:,i]=y1
            #np.savetxt('output/Sep28_vhhAcc_pyCAMBunlensed(fix)_scalar_dCls_'+param+'.csv',np.nan_to_num(python),delimiter=',')
            #np.savetxt('output/Sep28_vhhAcc_unfixKT_scalar_dCls_'+param+'.csv',np.nan_to_num(axion),delimiter=',')
        ax.plot(x1,y1,label='axion')
        ax.plot(x2,y2,label='python')
        ax.set_xlim([0,100])
        print paramNames[index]
        print y1[:7500]
        print y2[:7500]
        ax.set_title(param)
        plt.legend()
        #plt.show()
        #plt.savefig('pyCAMB_axion_'+powers[i]+'(vhhAccSep28).png',format='png')
    sys.exit()
    plt.close(fig)
'''
fig = plt.figure(figsize=(16,10))
data1 = np.loadtxt('output/Sep28_vhhAcc_pyCAMBunlensed_scalar_fCls.csv',delimiter=',')
#data1 = np.loadtxt('output/June7_newoptimal_vhhAcc_unlensed_scalar_fCls.csv',delimiter=',')
data2 = np.loadtxt('output/Sep28_vhhAcc_unfixKT_unlensed_axion_fCls.csv',delimiter=',')
data3 = np.loadtxt('output/Sep28_vhhAcc_fixKT_unlensed_axion_fCls.csv',delimiter=',')
powers = ['TT','EE','BB','TE','KK','KT'] 
for i in range(6):
    ax = fig.add_subplot(3,2,i+1)
    if i > 10:
        print 'wrong!'
        x1 = np.array(range(len(data1[:,i])))
        y1 = data1[:,i] * (x1*(x1+1.))
        x2 = np.array(range(len(data2[:,i])))
        y2 = data2[:,i] * (x2*(x2+1.))
        x3 = np.array(range(len(data3[:,i])))
        y3 = data3[:,i] * (x3*(x3+1.))
    else: 
        x1 = np.array(range(len(data1[:,i])))
        y1 = data1[:,i]
        x2 = np.array(range(len(data2[:,i])))
        y2 = data2[:,i]
        x3 = np.array(range(len(data3[:,i])))
        y3 = data3[:,i]
        if i == 5:
            y2 = y2 / np.sqrt(x2*(x2+1.)) * (2.7255*1.e6) * 2
            y1 = np.nan_to_num(data1[:,i] / np.sqrt(x1*(x1+1.)) *(2.7255*1.e6) * 2)
            #data3[:,i] = y3 * np.sqrt(x3*(x3+1.)) /(2.7255*1.e6) / 2
            data1[:,i] = y1

    #ax.plot(x1,y1,label='py - sep24 add YHe')
    #ax.plot(x1,y1,label='py - sep28')
    ax.plot(x2,y2,label='axion unfixKT vhh')
    ax.plot(x3,y3,label='axion fixKT vhh')
    ax.set_title(powers[i])
    #ax.set_xlim([0,1000])
    ax.set_xscale("log")
    #print y1[:7500]
    print y2[:7500]
    print y3[:7500]
    plt.legend(loc='upper right')
#np.savetxt('output/Sep28_vhhAcc_unfixKT_unlensed_axion_fCls.csv',data3,delimiter=',')
np.savetxt('output/Sep28_vhhAcc_pyCAMBunlensed(fix)_scalar_fCls.csv',data1,delimiter=',')
#plt.savefig('Sep28_pyCAMB_axion(fixKT_vhh).png',format='png')
plt.show()

