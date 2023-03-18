'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from . import partial_binodal, io
import time

# Function optimizer performs the optimization
def optimizer(mode,
              regime,
              dilute_array,
              dense_array,
              init_guess,
              free_parameter_checklist,
              temperature_offset,
              fitting_method,
              maxfev,
              xatol,
              fatol,
              partial_fit,                               
              silent=False):

    """
    Function for performing optimization of the phase diagram fit using the given theory.
    Input parameters are defined below. This function uses scipy.optimize.minimize for fitting.

    
    Parameter
    --------------------
    mode : str
        The name of the theory module to use. By default, fireball-fit and fireball-draw use flory_huggins_3B
        
    regime : str (either "low" or "high")
        If "low", then the algorithm assumes an upper critical system
        If "high", then the algorithm assumes a lower critical system
        By default, fireball-fit and fireball-draw use "low"
    
    dilute_array : np.ndarray (2D)
        Array of points along the dilute arm of the binodal. Includes
        both concentrations and temperature. Concentrations should be in
        volume fraction (i.e. between 0 and 1).
    
    dense_array : np.ndarray (2D)
        Array of points along the dense arm of the binodal. Includes
        both concentrations and temperature. Concentrations should be in
        volume fraction (i.e. between 0 and 1).

    init_guess : list
        List that defines initial guesses for the set of parameters for the given theory module. The default theory module, flory_huggins_3B, requires 
        3 parameters. Element 0 is the chi energy term (w1), element 1 is the three body interaction parameter (w3) and element 2 
        is the degree of polymerization (n).

    free_parameter_checklist : list
        List that defines which parameters should be fixed and which should be optimized. For example, if the default theory module, 
        flory_huggins_3B, is used, then a free_parameter_list of [1, 1, 0] suggests that w1 and w3 are optimized but n is kept fixed.

    temperature_offset : float
        Converts the y-axis units when plotting the data. Should be set to 273.15 to convert the units
        from Kelvin to degrees Celsius. Should bet set to 0 to keep the units in Kelvin.

    fitting_method : str
        Numerical fitting method that scipy.optimize will use. Default
        recommended = 'Nelder-Mead'.

    maxfev : int
        Max number of calls to the function (see scipy docs)

    xatol : float
        Absolute error in f(x) to meet for convergence in terms
        of fitted parameters (i.e. x tolerance)

    fatol : floal
        Absolute error in f(x) to meet convergence in terms of 
        residual output (i.e. f tolerance).

    partial_fit : bool
       Flag which, if set to true, will pop up a fit of the 
       parameters based on the partial fit

    silent : bool
        Flag which, if set to true, means the function prints
        nothing to STDOUT nor generates a plot (i.e. this overrides
        the partial_fit flag.


    Returns
    -----------------

    scipy.optimize.minimize object
    
        
    """

    # Determine time at start of fitting 
    start_time = time.time()              
    

    # Define and perform our fit
    # The parameters worth noting are "method," which by default is Nelder-Mead (a gradient-free method with fast convergence) and the options, which are specific to the method 
    if not silent:
        print('Starting fitting...')
    
    free_parameter_list = []
    fixed_parameter_list = []
    for index, param in enumerate(init_guess):
        if free_parameter_checklist[index] == 1:
            free_parameter_list.append(param)
        else:
            fixed_parameter_list.append(param)
            
    fit = minimize(partial_binodal.residual_finder, free_parameter_list, (fixed_parameter_list, free_parameter_checklist, mode, regime, dilute_array, dense_array, temperature_offset), method = fitting_method, options={'maxfev': maxfev, 'xatol': xatol, 'fatol': fatol})

    # If silent is true then do nothing
    if silent:
        pass
    else:

        # Print out fitting time and a summary of the fit
        delta = str(time.time() - start_time)
        print("\n\nFitting time: {}".format(delta))
        print('\n=======================================')
        print('Summary of fitting procedure:')
        print(fit)
        print('\n=======================================')

        # This line plots our data along with the last partial binodal produced during the fitting process in order to make sure that the residuals were calculated correctly
        # The "1" in the argument list instructs the function to show the plot at the end
        # If any of the points in the partial binodal appear incorrect, the parameters in config.py or in partial_binodal.py should be altered to gain more accuracy
        # This will also print the error list that describes how accurately the program is computing the binodal. In my experience, if the errors are ~10^-14 or less, the program is functioning correctly
        # If the errors are closer to 10^-10 - 10^-5, the parameters for computing the binodal may need to be modified
        if partial_fit:
            partial_binodal.residual_finder(fit.x, fixed_parameter_list, free_parameter_checklist, mode, regime, dilute_array, dense_array, temperature_offset, plotter=True, error=True, title='Partial fit')

    return fit
