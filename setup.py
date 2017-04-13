from distutils.core import setup, Extension
import os


setup(name='pyfisher',
      version='0.1',
      description='Modules for CMB Fisher Forecasts',
      author='Mathew Madhavacheril, Ho Nam Nguyen',
      packages=['pyfisher'],
      install_requires=['numpy',
                        'matplotlib',
                        'scipy'],
      zip_safe=False)
