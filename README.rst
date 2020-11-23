========
pyfisher
========


.. image:: https://img.shields.io/pypi/v/pyfisher.svg
        :target: https://pypi.python.org/pypi/pyfisher

.. image:: https://img.shields.io/travis/msyriac/pyfisher.svg
        :target: https://travis-ci.com/msyriac/pyfisher

.. image:: https://readthedocs.org/projects/pyfisher/badge/?version=latest
        :target: https://pyfisher.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


.. image:: https://pyup.io/repos/github/msyriac/pyfisher/shield.svg
     :target: https://pyup.io/repos/github/msyriac/pyfisher/
     :alt: Updates



Fisher forecasting for cosmological surveys

pyfisher is a python package for calculating Fisher matrices and for forecasting parameter uncertainties for cosmological surveys.

![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)  **Version 2 is a total revamp, so if you're used to what this software looked like before November 2020, you should switch to the `legacy` branch.** 

While the new version does not (yet) provide an interface for CMB lensing noise curves with iterative
delensing like the old one did, it has a simplified API, lots of pre-calculated Fishers, and a tool to reparametrize into a <img src="https://render.githubusercontent.com/render/math?math=\sigma_8">  parameterization.


* Free software: BSD license
* Documentation: https://pyfisher.readthedocs.io.


Installation
------------

Install in two steps:

1. Git clone the repository
2. `cd` into the repository and run `pip install -e . --user`.

The latter step just copies symbolic links to the relevant modules into a directory (managed by pip) that is in your python path.

Once this is done, you should be able to do things like

``
import pyfisher
``

from anywhere on your system.


Basic Usage
-----------

See and run `python bin/test_lensing.py planck` to reproduce Planck constraints and get a feel for how to use this package.

