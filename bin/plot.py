from __future__ import print_function
from chainconsumer import ChainConsumer



mean = [1, 5,6]
cov = [[1, 1, 0], 
       [1, 3, 0],
       [0, 0, 4]]
parameters = ["a", "b","c"]

c = ChainConsumer()
c.add_covariance(mean, cov, parameters=parameters, name="Cov")
c.add_marker(mean, parameters=parameters, name="Marker!", marker_style="*", marker_size=100, color="r")
c.configure(usetex=False, serif=False,sigma2d=True,sigmas=[1])
fig = c.plotter.plot()

fig.set_size_inches(3 + fig.get_size_inches()) 
fig.savefig('test.png')
