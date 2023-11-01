FIREBALL (FItting REal BinodALs with non-Linear models)
==============================

## About
A general-purpose Python package for fitting and building phase diagrams. Also includes methods for fitting coil-to-globule transitions.

## Installation
As of right now FIREBALL can be installed from the `github` repository. Specifically, this involves;

1. Cloning the repo locally - i.e.

		git clone https://github.com/Papplab/FIREBALL.git
		
2. Move into the main FIREBALL directory  (where `setup.py` is) and run

		pip install .
		
This will install FIREBALL locally. Alternatively, you can run

	pip install -e .
  
This will install FIREBALL locally but link to this specific location of the code, so changes made are instantaneously reflected in the global install

## Command-line usage
Fireball comes with multiple command-line tools. The two main tools used for plotting binodal data are `fireball-fit` and `fireball-draw`. There are also two tools used for plotting single chain coil-to-globule transitions, which are `fireball-single-chain-fit` and `fireball-single-chain-draw`

### fireball-fit
`fireball-fit` is the main command for fitting to binodals. The simplest usage is

	fireball-fit --filename <name of csv or .xls file with dense/dilute phase data>
	
> **Note** `pandas` deprecated the ability to read `.xlsx` files late in 2020, hence the `.xls` requirement).

This will fit the provided dense/dilute phase using the default options (discussed below).

#### Input files
Input files should follow the following format:

	Temp, Dense, Dilute, Dense_error, Dilute_Error
	
Note that a header line matching this is required, and then each subsequent line should populate `Temp` and at least one of `Dense` or `Dilute`. Concentrations should be in mg/mL

For example, a valid `.csv` file might look like:

	Temp, Dense, Dilute, Dense_Error, Dilute_Error
	277.15, , 0.169, , 0.003
	278.15, 393, 0.177, , 0.004
	279.15, , 0.207, , 0.004
	281.15, , 0.304, , 0.004

#### Input options
`fireball-fit` offers the following command line options:
	    
	  -h, --help            show this help message and exit
	  --filename FILENAME   Input file with data for fitting
	  --initial_guess INITIAL_GUESS
	                        [OPTIONAL] Initial guess for values (format X_Y_Z
	                        where X,Y,Z depend on the fitting procedure.
	                        [Default = 209_0.31_100]
	  --free_param_list FREE_PARAMETER_LIST
	                        [OPTIONAL] Define which parameters should be free vs
	                        fixed (format X_Y_Z where X,Y,Z are 0 if the parameter
	                        is fixed and 1 if it is free. [Default = 1_1_0]
	  --conversion-factor CONVERSION_FACTOR
	                        [OPTIONAL] Conversion factor that converts
	                        concentration in mg/mL to volume fraction [Default =
	                        1310 mg/mL]
	  --temperature-offset TEMPERATURE_OFFSET
	                        [OPTIONAL] Value that defines the temperature offset
	                        for plotting [Default = 273.15, assumes data should
	                        be plotted in Celsius]
	  --LCST-check LCST_CHECK
	                        [OPTIONAL] Flag which, if set to true, tells the
	                        algorithm that this is a lower critical system [Default
	                        = False]
	  --fitting-method FITTING_METHOD
	                        [OPTIONAL] Method used to fit [Default = 'Nelder-Mead']
	  --maxfev MAXFEV       [OPTIONAL] Max number of function evaluations allowed
	                        during optimization [Default=200]
	  --xatol XATOL         [OPTIONAL] The absolute error in the parameter set
	                        between iterations that is acceptable for convergence
	                        [Default=0.1]
	  --fatol FATOL         [OPTIONAL] The absolute error in the square-residuals
	                        between iterations that is acceptable for convergence
	                        [Default=0.1]
	  --increment INCREMENT [OPTIONAL] The phase diagram increment [Default=1.0]
	  --mode MODE           [OPTIONAL] Defines the fitting mode to be used.
	  					  Valid options are: 'flory-huggins-3B' and 'GCT'
	  					  [Default=flory-huggins-3B]
	  --outname OUTNAME     [OPTIONAL] Prefix of output files where
	  					  binodal/spinodal data are written [Default = Fit]
	  --outname-fitfile OUTNAME_FITFILE
	                        [OPTIONAL] Name of output file where final fitting
	                        parameters are written
	                        [Default = fit_params.csv]
	  --show-partial-fit    Flag which, if set to true, will show the binodal
	                        generated from the partial fit.
	  --noplot              Flag which, if provided, means no figure is generated
	                        to screen but the binodal data are just saved to disk
	  --version             Flag which, if set to true, prints software version
	  					  and exits the program

A few non-intuitive things:

* Right now the only `--mode` options are `flory-huggins-3B`, which uses a modified version of Flory-Huggins theory, and `GCT`, which uses Gaussian Cluster theory. The Flory-Huggins theory is the same as that used in Martin et al. Science (2020). The Gaussian Cluster theory is based on Raos and Allegra. J Chem Phys (1996) and is well-described in Zeng et al. Biophysical Journal (2020).

* By default, the first two parameters are free and third parameter is fixed. These correspond to w1, w3, and n, respectively. In general, the value of a fixed parameter is kept the same as whatever is passed to `--inital_guess`. For example, if `--inital_guess` is set to `200_0.01_150` and `--free_param_list` is set to `0_1_1` then w1 will be fixed at 200, whereas w3 and n will be free parameters but initialized at 0.1 and 150, respectively.

* `--fitting-method` has only been tested with `Nelder-Mead`, so changing this flag may result in unexpected behavior.

* `--version` prints the FIREBALL version, which means `--version` from all scripts will give the same version. Versioning is done automatically by git.

As output this generates a few files. Assuming `--noplot` is not set, the function by default generates a plot to the screen and also saves the binodals in normal and semi-logspace (wrt phi) as PDFs. It will also save three `.csv` files, which by default are called `fit_binodal.csv`, `fit_spionodal.csv`, and `fit_params.csv`. These contain the binodal information, the spinodal information, and the summary of the fitted parameters, respectively.

### fireball-draw
`fireball-draw` lets you build and plot the data for a phase diagram by passing in parameters directly (rather than fitting to data). 

#### Input options
`fireball-draw` offers the following command line options:

	  -h, --help            show this help message and exit
	  --filename FILENAME   [OPTIONAL] Input file with information for plotting data
	  --param_list PARAM_LIST
	                        [OPTIONAL] Parameter values (format X_Y_Z where X,Y,Z
	                        depend on the fitting procedure. [Default =
	                        209_0.31_100]
	  --conversion-factor CONVERSION_FACTOR
	                        [OPTIONAL] Conversion factor that converts
	                        concentration in mg/mL to volume fraction [Default =
	                        1310 mg/mL]
	  --temperature-offset TEMPERATURE_OFFSET
	                        [OPTIONAL] Value that defines the temperature offset
	                        for plotting [Default = 273.15, assumes data should
	                        be plotted in Celsius]
	  --LCST-check LCST_CHECK
	                        [OPTIONAL] Flag which, if set to true, tells the
	                        algorithm that this is a lower critical system [Default
	                        = False]
	  --increment INCREMENT [OPTIONAL] The phase diagram increment [Default=1.0]
	  --mode MODE           [OPTIONAL] Defines the fitting mode to be used.
	  					  Valid options are: 'flory-huggins-3B' and 'GCT'
	  					  [Default=flory-huggins-3B]
	  --outname OUTNAME     [OPTIONAL] Prefix of output files where
	  					  binodal/spinodal data are written [Default = Draw]
	  --noplot              Flag which, if provided, means no figure is generated
	                        to screen but the binodal data are just saved to disk
	  --version             Flag which, if set to true, prints software version
	  					  and exits the program

Note that the `--filename` is optional, and the simplest usage

	fireball-draw
	
should generate a binodal to the screen.

### fireball-single-chain-fit
`fireball-single-chain-fit` is the main command for fitting to coil-to-globule transitions. It functions similarly to `fireball-fit` with a few important differences outlined below.

#### Input files
Input files should follow the following format:

	Temp, Rg, Rg_error
	
Note that a header line matching this is required, and then each subsequent line must populate `Temp` and `Rg`. 

For example, a valid `.csv` file might look like:

	Temp, Rg, Rg_Error
	277.15, 22.1, 0.4
	278.15, 24.1, 0.1
	279.15, 25.5, 0.2
	281.15, 25.8, 0.4

#### Input options
`fireball-single-chain-fit` offers the following command line options:


	  -h, --help            show this help message and exit
	  --filename FILENAME   Input file with data for fitting
	  --initial_guess INITIAL_GUESS
	                        [OPTIONAL] Initial guess for values (format X_Y_Z
	                        where X,Y,Z depend on the fitting procedure.
	                        [Default = 345_.02_.0076_3000]
	  --free_param_list FREE_PARAMETER_LIST
	                        [OPTIONAL] Define which parameters should be free vs
	                        fixed (format X_Y_Z where X,Y,Z are 0 if the parameter
	                        is fixed and 1 if it is free. [Default = 1_1_1_1]
	  --temperature-offset TEMPERATURE_OFFSET
	                        [OPTIONAL] Value that defines the temperature offset
	                        for plotting [Default = 273.15, assumes data should
	                        be plotted in Celsius]
	  --LCST-check LCST_CHECK
	                        [OPTIONAL] Flag which, if set to true, tells the
	                        algorithm that this is a lower critical system [Default
	                        = False]
	  --fitting-method FITTING_METHOD
	                        [OPTIONAL] Method used to fit [Default = 'Nelder-Mead']
	  --maxfev MAXFEV       [OPTIONAL] Max number of function evaluations allowed
	                        during optimization [Default=200]
	  --xatol XATOL         [OPTIONAL] The absolute error in the parameter set
	                        between iterations that is acceptable for convergence
	                        [Default=0.1]
	  --fatol FATOL         [OPTIONAL] The absolute error in the square-residuals
	                        between iterations that is acceptable for convergence
	                        [Default=0.1]
	  --increment INCREMENT [OPTIONAL] The temperature increment [Default=1.0]
	  --mode MODE           [OPTIONAL] Defines the fitting mode to be used.
	  					  Valid options are: 'GCT_Single_Chain'
	  					  [Default=GCT_Single_Chain]
	  --outname OUTNAME     [OPTIONAL] Name of output file where data are written
	  					  [Default = single-chain.csv]
	  --outname-fitfile OUTNAME_FITFILE
	                        [OPTIONAL] Name of output file where final fitting
	                        parameters are written
	                        [Default = fit_params.csv]
	  --show-partial-fit    Flag which, if set to true, will show the coil-
	  					  globule transition generated from the partial fit
	  --noplot              Flag which, if provided, means no figure is generated
	                        to screen but the single chain data are just saved
	                        to disk
	  --version             Flag which, if set to true, prints software version
	  					  and exits the program

As output this generates a few files. Assuming `--noplot` is not set, the function by default generates a plot to the screen and also saves the single chain data as a PDF. It will also save two `.csv` files, which by default are called `single-chain.csv` and `fit_params.csv`. These contain the single chain data and the summary of the fitted parameters, respectively.

### fireball-single-chain-draw
`fireball-single-chain-draw` lets you build and plot the data for a single chain coil-to-globule transition by passing in parameters directly (rather than fitting to data). It functions similarly to `fireball-draw` with a few important differences outlined below.

#### Input options
`fireball-single-chain-draw` offers the following command line options:


	  -h, --help            show this help message and exit
	  --filename FILENAME   [OPTIONAL] Input file with information for plotting data
	  --param_list PARAM_LIST
	                        [OPTIONAL] Parameter values (format X_Y_Z where X,Y,Z
	                        depend on the fitting procedure. [Default =
	                        345_.02_.0076_3000]
	  --temperature-offset TEMPERATURE_OFFSET
	                        [OPTIONAL] Value that defines the temperature offset
	                        for plotting [Default = 273.15, assumes data should
	                        be plotted in Celsius]
	  --LCST-check LCST_CHECK
	                        [OPTIONAL] Flag which, if set to true, tells the
	                        algorithm that this is a lower critical system [Default
	                        = False]
	  --increment INCREMENT [OPTIONAL] The temperature increment [Default=1.0]
	  --mode MODE           [OPTIONAL] Defines the fitting mode to be used.
	  					  Valid options are: 'GCT_Single_Chain'
	  					  [Default=GCT_Single_Chain]
	  --outname OUTNAME     [OPTIONAL] Name of output file where
	  					  single chain data are written [Default = 'single-chain.csv]
	  --noplot              Flag which, if provided, means no figure is generated
	                        to screen but the single chain data are just saved to disk
	  --version             Flag which, if set to true, prints software version
	  					  and exits the program

## Package usage
In addition, FIREBALL also lets you work directly inside Python. For example:

	from matplotlib import pyplot as plt
	from matplotlib.pyplot import figure
	from fireball import full_binodal
	
	# define some params
	W1=200
	W3=0.1
	n = 200
	mode = 'flory_huggins_3B'
	regime = 'low'
	
	# build figure frame
	figure(num=None, figsize=(5, 2), dpi=150, facecolor='w', edgecolor='k')
	ax = plt.gca()
	
	for x in range(150,250,20):
	    A = full_binodal.build_theoretical_binodal_flory_huggins_3B([x,  W4, n], mode, regime)
	    plt.plot(A[1],A[2],'k--', linewidth=0.5)
	    plt.plot(A[1],A[2],'.', alpha=0.5, markersize=2, label='%i'%(x))
	
	plt.xlabel(r'$\phi$')
	plt.ylabel(r'T')
	plt.legend(frameon=False)

## Architecture and development
General structure of the code is as follows:

All source code can be found in 

	fireball/

All theoretical functions, such as Flory-Huggins free energy formulations can be found in 

	fireball/theory

All configuration settings for FIREBALL can be found in 

	fireball/config.py

In addition to optimization hyper-parameters, this module also tells FIREBALL the starting temperatures for calculating full binodals and how to label various plots. Depending on your system of interest, you may need to change some of the values in this module.

All command-line tools can be found in

	scripts/
	
This is where I recommend starting to read through the code.

### Tutorial
We have included two short tutorials in this package to illustrate the process of fitting phase diagrams (tutorial 1) and coil-to-globule transitions (tutorial 2) in the `tutorials/` directory. Here, we will run through the process of fitting binodal data as an example to show how FIREBALL functions. After parsing input options the `fireball-fit` program does the following:

#### Stage 1 - fit the data
It calls `optimizer.optimizer()`, which is the wrapper code for fitting binodal data. 

`optimizer.optimizer()` is inside `fireball/optimizer.py`, which itself optimizes by applying `scipy.optimize.minimize` on the function `partial_binodal.residual_finder()`. 

`partial_binodal.residual_finder()` is inside `fireball/partial_binodal.py` and takes as input a set of energy parameters and datapoints, calculates a theoretical binodal using the given parameters, and returns a rescaled sum of residuals between the calculated binodal and the datapoints to the `scipy.optimize.minimize` function. In so doing, this function calls `partial_binodal_loop.partial_binodal_loop()`. For now, the residuals are calculated using a log-scale, to take into account the fact that protein phase diagram concentrations typically change by orders of magnitude, but this can be changed if needed.

`partial_binodal_loop.partial_binodal_loop()` is inside `fireball/partial_binodal_loop.py` and does the heavy lifting to calculate the theoretical binodal. This function loads components of the actual physical model from `fireball/theory/flory_huggins_3B.py` (the default theory module).

If performed in the command-line, `fireball-fit` will print out the rescaled sum of residuals after each loop through `partial_binodal.residual_finder()`. Eventually, the `scipy.optimize.minimize` function will either converge on an optimized set of parameters or perform the maximum number of function evaluations.

#### Stage 2 - plot the data
`fireball-fit` then calls `full_binodal.binodal_maker()`, which is inside `fireball/full_binodal.py`, using the newly fitted parameters, which plots the full binodal alongside the initial data set. The binodal, spinodal, and fitted parameters will also be saved to disk

### Adding new theory modules
New theories can be readily added to FIREBALL by creating a new file in the `fireball/theory` directory. I recommend using an existing theory as a template. If the theory has a similar overall structure as Flory-Huggins theory, then you should only need to add the name of the new theory module to the list of allowed theories in the script files, and FIREBALL should just work. If the theory has a unique overall structure, then it may require changes to the underlying source code.

## Troubleshooting
Depending on the data and the initial parameters, it is possible that FIREBALL will not adequately calculate phase diagrams. Hopefully, Python spits out a reasonable error message that informs you of the issue. In general, you can take the following steps to try to troubleshoot the issue:

### Step 1: Change initial parameters
First, change the initial parameters. Often, the initial parameters are either unphysical or too far from the optimized parameters for FIREBALL to function correctly. I recommend using `fireball-draw` to determine how the binodal calculated from your initial parameters compares to the data set in the absence of any fitting.

### Step 2: Change the `config.py` module
Second, if changing the initial parameters did not help, try modulating values in the `config.py` module. In particular, make sure that your starting binodal temperature is not too high or too low. The first five variable in this module may also be of interest, in particular "partial_binodal_increment", which may need to be decreased.

### Step 3: Change the bounds/initial guesses when calculating binodal values in `fireball/partial_binodal_loop.py`
If the above steps did not help, it is possible that you will need to change the source code. Most likely, the algorithm is having trouble calculating correct binodal values in the `fireball/partial_binodal_loop.py` module. If this seems to be the case, I recommend playing with the initial guess and bound variables in this module. These are named "x0" and "bounds", respectively. Unfortunately, you will likely need to familiarize yourself with the code and understand the issue at hand to adequately change these variables. Depending on the issue, you may need to make similar changes in the `fireball/full_binodal.py` module to correctly calculate the full binodal at the end of the script.

## Development

### Copyright

Copyright (c) 2019-2023, Pappu lab

#### Acknowledgements
 
The FIREBALL code was written by Mina Farag and Alex Holehouse
