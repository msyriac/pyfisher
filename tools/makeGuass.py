import numpy as np
import matplotlib.pyplot as plt


def gaussian(x, mu, sig):
    return np.exp(-(x-mu)**2/(2*sig**2))/np.sqrt(2*np.pi*sig**2)

mu = 1.4
sig = 0.2
N=2000

fig = plt.figure()
ax = fig.add_subplot(111)
x = np.linspace(0,3.5,N)
y = gaussian(x,mu,sig)

data = np.hstack([x.reshape(N,1),y.reshape(N,1)])
print x[1]-x[0],y
print (y*(x[1]-x[0])).sum()
print x,y,data

ax.plot(x,y)
#ax.set_xlim([60,80])
plt.show()
np.savetxt('output/dndz_Das_lensing.csv',data,delimiter=',')
