# Forecasting Tutorial

The recommended procedure is to use

```
python bin/driver.py input/external.ini
```

for external Fisher matrices (Planck, BAO, etc.) so that those get saved to `output/`.

Then use

```
python bin/lensing.py input/params.ini
```

to create a Fisher matrix for the primary CMB and for Clkk and save this to `output/`.

Then you can use whatever combinations of the saved Fisher matrices in a separate script.

## External Fisher matrices

Typically, you will want to do something like this:

1. Include Planck TT only for l<30 over fsky=.65
2. Include Planck TT,TE,EE for 30<l<2500 over fsky=0.65 - fsky_exp, where fksy_exp is the fsky of the experiment you are forecasting for
3. Include a BOSS or DESI BAO Fisher matrix
4. Include a tau prior of 0.01 from Planck

You are encouraged to use bin/driver.py with input/external.ini to achieve 1-3. Manually add a tau prior yourself to the sum of your Fisher matrices at the end. driver.py lets you save each included section of input/external.ini as a txt file whose name you specify in external.ini. Please edit the commented lines in input/external.ini appropriately.