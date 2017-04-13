import numpy as np

row = 160
column = 8
root = '/gpfs/scratch/nhnguyen/testStepSizeDump/'
#testList = ['H0','ombh2','omch2','tau','ns','As','mnu','w']
testList = ['mnu']
for testParam in testList:
    name = 'June23_testStepSize_240meV_unlensed_scalar_'+testParam+'_LCDM+mnu_S4'
    data = np.zeros([row,column])
    for i in range(row):
        data[i,:] = np.loadtxt(root+name+'_part'+str(i+1)+'.csv',delimiter=',')
    #print data
    np.savetxt(root+name+'.csv',data,delimiter=',')
    print "Saved file",root+name+'.csv'
'''
row = 5
column = 8
root = '/gpfs/scratch/nhnguyen/testStepSizeDump/'
data = []
for i in range(row):
    newdata = np.loadtxt(root+'June4_'+str(i+1)+'_testStepSize_veryhighhighAcc_unlensed_scalar_mnu_LCDM+mnu_optimal_CMB+BAO.csv',delimiter=',')
    #print len(data)
    if len(data) == 0:
        data = newdata
    else:
        data = np.concatenate((data,newdata))
    #print data
#print len(data),data[:,0]
np.savetxt(root+'June4_testStepSize_veryhighhighAcc_unlensed_scalar_mnu_LCDM+mnu_optimal_CMB+BAO.csv',data,delimiter=',')

'''
