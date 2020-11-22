# pyfisher 2

pyfisher is a python package for calculating Fisher matrices and for forecasting parameter uncertainties for cosmological surveys.

```diff
- Version 2 is a total revamp, so if you're used to what this software looked like before November 2020, you should
switch to the `legacy` branch. 
```

While the new version does not (yet) provide an interface for CMB lensing noise curves with iterative
delensing like the old one did, it has a simplified API, lots of pre-calculated Fishers, and a tool to reparametrize into a <img src="https://render.githubusercontent.com/render/math?math=\sigma_8">  parameterization.


## Installation

You will need to install a dependency called `orphics` that is not available on pip, conda, etc.
https://github.com/msyriac/orphics

The installation procedure for that dependency is the same as for this package. Repeat the following for each of `orphics` and `pyfisher` separately.

1. Git clone the repository
2. `cd` into the repository and run `pip install -e . --user`.

The latter step just copies symbolic links to the relevant modules into a directory (managed by pip) that is in your python path.

Once this is done, you should be able to do things like

``
import pyfisher
``

from anywhere on your system.


## Basic Usage

See and run `python bin/test_lensing.py planck` to reproduce Planck constraints and get a feel for how to use this package.
