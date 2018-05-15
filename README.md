# PyFisher

PyFisher is a python package for calculating Fisher matrices and forecasting parameter uncertainties for CMB experiments. 

## Installation

You will need to install a dependency called `orphics` that is not available on pip, conda, etc.
https://github.com/msyriac/orphics

The installation procedure for that dependency is the same as for this package. Repeat the following for each of `orphics` and `pyfisher` separately.

1. Git clone the repository
2. `cd` into the repository and run `pip install -e . --user`.

The latter step just copies symbolic links to the relevant modules into a directory (managed by pip) that is in your python path.

Once this is done, you should be able to do things like

``
import pyfisher.clFisher as clFish
``

from anywhere on your system.


## Basic Usage

1. Edit `input/params.ini`. It is documented. You can include or exclude BAO by changing `otherFishers`, and the details of lensing reconstruction by changing the `lensing` section. There are a whole bunch of sections for different experiment configurations too.
2. See how `tests/testFisher.py` uses this ini file. Inside the root repo directory, run
``
python tests/testFisher.py <experiment section name> <lensing section name>
``
For example
``
python tests/testFisher.py S4-5m lensing
``
which will tell you what the mnu constraint and lensing S/N is. You can adapt this script to, for example, override one of the ini configurations and make plots of constraints and S/N varying that configuration.

## Advanced Usage


### Make derivatives

Edit `input/makeDefaults.ini` and specify a fiducial cosmology, derivative step sizes and a prefix under which to save the derivatives to  `output/`. Run with

`python bin/makeDerivs.py`

### Re-make Planck or BAO Fisher matrices.

Edit `input/fisherDebug.ini` and run
```
python bin/driver.py input/fisherDebug.ini
```

### A tutorial

A very common forecasting scenario is combining a future CMB experiment with existing Planck data and BAO. Please read this tutorial for suggestions on how to do this with pyfisher.

https://github.com/msyriac/pyfisher/blob/master/forecasting.md
