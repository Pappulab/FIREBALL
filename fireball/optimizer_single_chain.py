'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from . import partial_single_chain
import time

# Function optimizer performs the optimization
def optimizer(mode,
              data_array,
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
        The name of the theory module to use. By default, fireball-single-chain-fit and fireball-single-chain-draw use GCT_Single_Chain
        
    data_array : np.ndarray (2D)
        Array of points along the coil-globule transition. Includes both Rg and temperature.
    
    init_guess : list
        List that defines initial guesses for the set of parameters defined for the given theory module.
        The default theory module is GCT_Single_Chain.

    free_parameter_checklist : list
        List that defines which parameters should be fixed and which should be optimized.

    temperature_offset : float
        Converts the y-axis units when plotting the data. Should be set to 273.15 to convert the units
        from Kelvin to degrees Celsius. Should bet set to 0 to keep the units in Kelvin.

    fitting_method : str
        Numerical fitting method that scipy.optimize will use. Default recommended = 'Nelder-Mead'.

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
    # The parameters worth noting are "method," which by default is Nelder-Mead (a gradient-free method
    # with fast convergence) and the options, which are specific to the method 
    if not silent:
        print('Starting fitting...\n')
    
    free_parameter_list = []
    fixed_parameter_list = []
    bounds_list = []
    for index, param in enumerate(init_guess):
        if free_parameter_checklist[index] == 1:
            free_parameter_list.append(param)
            if index == 2:
                # Place a lower bound on the three-body interaction parameter
                bounds_list.append((10 ** -10, None))
            else:
                bounds_list.append((None, None))
        else:
            fixed_parameter_list.append(param)
    
    # Perform the fit
    fit = minimize(partial_single_chain.residual_finder, free_parameter_list,
                   (fixed_parameter_list, free_parameter_checklist, mode, data_array,
                    temperature_offset), method = fitting_method, bounds = bounds_list,
                   options={'maxfev': maxfev, 'xatol': xatol, 'fatol': fatol})

    if fitting_method == 'Nelder-Mead':
        print('\nEstimating covariance matrix...')
        
        final_simplex = fit.final_simplex
        final_params = final_simplex[0]
        final_funcs = final_simplex[1]
        final_centroid = np.mean(final_params, axis = 0)
        final_centroid_func = partial_single_chain.residual_finder(final_centroid, fixed_parameter_list, free_parameter_checklist, mode, data_array, temperature_offset, silent=True)
        threshold = final_centroid_func + final_centroid_func * .01
        new_params = final_simplex[0]
        num_params = len(free_parameter_list)
        num_vertices = num_params + 1
        new_funcs = np.zeros(num_params + 1)
        
        threshold_check = True
        while threshold_check:
            threshold_check = False
            new_params = new_params + 2 * (new_params - final_centroid)
            for index, param_set in enumerate(new_params):
                new_funcs[index] = partial_single_chain.residual_finder(param_set, fixed_parameter_list, free_parameter_checklist, mode, data_array, temperature_offset, silent=True)
            for new_func in new_funcs:
                if new_func < threshold:
                    threshold_check = True
        
        new_params_array = [[] for i in range(num_vertices)]
        new_func_matrix = np.zeros((num_vertices, num_vertices))
        for i in range(num_vertices):
            for j in range(i + 1):
                new_params_array[i].append(np.mean([new_params[i], new_params[j]], axis = 0))
                new_func_matrix[i, j] = partial_single_chain.residual_finder(new_params_array[i][j], fixed_parameter_list, free_parameter_checklist, mode, data_array, temperature_offset, silent=True)
        
        a0 = new_func_matrix[0, 0]
        a_vector = np.zeros(num_params)
        for i in range(num_params):
            a_vector[i] = 2 * new_func_matrix[i + 1, 0] - 0.5 * (new_func_matrix[i + 1, i + 1] + 3 * new_func_matrix[0, 0])
        
        b_matrix = np.zeros((num_params, num_params))
        for i in range(num_params):
            for j in range(num_params):
                if j > i:
                    b_matrix[i, j] = 2 * (new_func_matrix[j + 1, i + 1] + new_func_matrix[0, 0] - new_func_matrix[i + 1, 0] - new_func_matrix[j + 1, 0])
                else:
                    b_matrix[i, j] = 2 * (new_func_matrix[i + 1, j + 1] + new_func_matrix[0, 0] - new_func_matrix[i + 1, 0] - new_func_matrix[j + 1, 0])
        b_inv_matrix = np.linalg.inv(b_matrix)
        
        q_matrix = np.zeros((num_params, num_params))
        for i in range(num_params):
            for j in range(num_params):
                q_matrix[j, i] = new_params[i + 1][j] - new_params[0][j]
        
        c_matrix = np.matmul(np.matmul(q_matrix, b_inv_matrix), np.transpose(q_matrix))
        
        y_min = a0 - np.matmul(np.matmul(a_vector, b_inv_matrix), np.transpose(a_vector))
        num_points = data_array.shape[0]
        sigma2 = y_min / (num_points - num_params)
        c_matrix_final = c_matrix * 2 * sigma2

    # If silent is true then do nothing
    if not silent:

        # Print out fitting time and a summary of the fit
        delta = str(time.time() - start_time)
        print("\n\nFitting time: {}".format(delta))
        print('\n=======================================')
        print('Summary of fitting procedure:')
        print(fit)
        print('\nEstimated covariance matrix:')
        print(c_matrix_final)
        print('=======================================')

        # This line plots our data along with the last partial coil-globule transition produced
        # during the fitting process in order to make sure that the residuals were calculated correctly.
        # The "1" in the argument list instructs the function to show the plot at the end.
        # If any of the points appear incorrect, the parameters in config.py or in partial_binodal.py
        # should be altered to gain more accuracy.
        # This will also print the error list that describes how accurately the program is computing
        # the coil-globule transition.
        # If the errors are large, the parameters for computing the coil-globule transition may need to be modified
        if partial_fit:
            partial_single_chain.residual_finder(fit.x, fixed_parameter_list, free_parameter_checklist, mode, data_array, temperature_offset, plotter=True, error=True, title='Partial fit')

    return fit
