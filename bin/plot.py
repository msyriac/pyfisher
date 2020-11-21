from __future__ import print_function
from chainconsumer import ChainConsumer


def contour_plot(fisher,fiducials,fname,name=None):
    from chainconsumer import ChainConsumer
    mean = [fiducials[key] for key in fisher.params]
    cov = fisher.values
    parameters = fisher.params

    c = ChainConsumer()
    c.add_covariance(mean, cov, parameters=parameters, name=name)
    c.add_marker(mean, parameters=parameters, marker_style="*", marker_size=100, color="r")
    c.configure(usetex=False, serif=False,sigma2d=True,sigmas=[1])
    fig = c.plotter.plot()
    fig.set_size_inches(3 + fig.get_size_inches()) 
    fig.savefig(fname)
