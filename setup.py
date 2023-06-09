"""
fireball
Python package for fitting FH data
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name='fireball',
    author='Pappu lab',
    author_email='alex.holehouse@wustl.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='LGPLv3',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    install_requires=[
        'matplotlib',
        'pandas',
        'numpy',
        'scipy',
        'seaborn',
        'PyYAML'
    ],

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # This provides an OS-agnostic means of installing executable scripts.
    entry_points = {
        'console_scripts': [
            'fireball-draw = fireball.scripts.fireball_draw:main',
            'fireball-fit = fireball.scripts.fireball_fit:main',
            'fireball-single-chain-draw = fireball.scripts.fireball_draw:main',
            'fireball-single-chain-fit = fireball.scripts.fireball_single_chain_fit:main',
        ],
    },

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    url='https://github.com/Pappulab/FIREBALL',  # Website

    platforms=['Linux',
               'Mac OS-X',
               'Unix',
               'Windows'],            # Valid platforms your code works on, adjust to your flavor

    python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    zip_safe=True,

)
