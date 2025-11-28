"""
fireball
Python package for fitting FH data
"""

import os

# Generate _version.py if missing and in the Read the Docs environment
if os.getenv("READTHEDOCS") == "True" and not os.path.isfile('../fireball/_version.py'):   
    import versioningit            
    __version__ = versioningit.get_version('../')
else:
    from ._version import __version__

# Add imports here
from fireball import *

