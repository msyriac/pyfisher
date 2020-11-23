"""Top-level package for pyfisher."""

__author__ = """Mathew Madhavacheril"""
__email__ = 'mathewsyriac@gmail.com'
__version__ = '2.0.0'

from .pyfisher import *

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
