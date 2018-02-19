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

