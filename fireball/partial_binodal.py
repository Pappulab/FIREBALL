'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

# import required packages
import importlib
import numpy as np
import pandas as pd
from scipy import optimize

from .fireballexceptions import FireballException
from .partial_binodal_loop import partial_binodal_loop
from .array_fixer import array_fixer
from . import config

def residual_finder(free_parameter_list, fixed_parameter_list, free_parameter_checklist,
                    mode, regime, dilute_array, dense_array, temperature_offset,
                    plotter=False, error=False, title=''):
    """

    This program is used to construct a partial binodal using a given set of parameters and a
    theory module. The default theory module is flory_huggins_3B.py, which involves a
    Flory-Huggins model with a 3-body correction term. Parameters can be designated as
    "fixed" or "free" (to be optimized) via the free_parameter_list argument (see below).
    
    This program finds the points along the binodal that correspond to the independent
    variable (temp, salt, etc.) values of the data. The residuals for these data points
    are the difference in phi values between the data points and the points on the
    partial binodal. This program can be used in conjunction with Optimizer.py in order
    to find the best fit parameters.

    This function is written such that it is compatible with scipy.optimize.minimize,
    i.e. first argument is a vector of parameters to be optimized over. This function is
    not really much use OTHER than in the context of being passed to minimize, which happens 
    in the optimizer.py module.

    Parameters
    ------------
    free_parameter_list : list
        List that defines the set of parameters set to be optimized for the given theory
        module. The default theory module, flory_huggins_3B, requires 3 parameters.
        Element 0 is the chi energy term (w1), element 1 is the three body interaction
        parameter (w3) and element 2 is the degree of polymerization (n).

    fixed_parameter_list : list
        List that defines the set of parameters set to be fixed for the given theory module.
        The default theory module, flory_huggins_3B, requires 3 parameters. Element 0 is
        the chi energy term (w1), element 1 is the three body interaction parameter (w3)
        and element 2 is the degree of polymerization (n).

    free_parameter_checklist : list
        List that defines which parameters should be fixed and which should be optimized.
        For example, if the default theory module, flory_huggins_3B, is used, then a
        free_parameter_list of [1, 1, 0] suggests that w1 and w3 are optimized but n
        is kept fixed.
    
    mode : str
        The name of the theory module to use. By default, fireball-fit and fireball-draw
        use flory_huggins_3B.
        
    regime : str (either "low" or "high")
        If "low", then the algorithm assumes an upper critical system
        If "high", then the algorithm assumes a lower critical system
        By default, fireball-fit and fireball-draw use "low"

    dilute_array : np.ndarray (3xm)
        A 2D array where column 0 is the dilute phase concentration, column 1
        is the temperature of the binodal, and column 3 is the error in the concentration.
        This is automatically generated as output
        element 0 from the fireball.io.read_file() function

    dense_array : np.ndarray (3xm)
        A 2D array where column 0 is the dense phase concentration, column 1
        is the temperature of the binodal, and column 3 is the error in the concentration.
        This is automatically generated as output
        element 1 from the fireball.io.read_file() function

    temperature_offset : float
        Converts the y-axis units when plotting the data. Should be set to 273.15 to convert the units
        from Kelvin to degrees Celsius. Should bet set to 0 to keep the units in Kelvin.

    plotter : bool
        Flag which, if set to true, means that the partial binodal fit is plotted to screen
        (only). Useful when debugging a fit, but in general not necessary. Default = False.

    error : bool
        Flag which, if set to true, prints out the list of errors in calculating the
        binodal points that correspond to given data points. This can be helpful in
        troubleshooting whether the algorithm is correctly computing binodal points.
        In my experience, error values should be in the 10^-20 - 10^-15 range.
        Values greater than 10^-10 may warrant a change in the parameters in config.py or 
        partial_binodal_loop.py that are used to create upper/lower bounds and guesses
        while determining the binodal
    
    title : str
        If plotter=True this lets you pass a title to the displayed plot

    Returns
    --------
    float
       The only thing this function returns is a float, which corresponds to the normalized
       residuals of the fit
    
    """
    
    # Import the correct theory module
    theory_module = importlib.import_module('fireball.theory.' + mode)
    
    # Define which parameters are fixed vs. free, as well as initial guesses for each
    parameter_list = []
    counter_free = 0
    counter_fixed = 0
    for element in free_parameter_checklist:
        if element == 1:
            parameter_list.append(free_parameter_list[counter_free])
            counter_free += 1
        else:
            parameter_list.append(fixed_parameter_list[counter_fixed])
            counter_fixed += 1
    print('Current parameter list is ' + str(parameter_list))
    
    # Calculate the critical point. Method varies depending on the theory used
    if mode == 'GCT':
        # The second derivative of the free energy should have two zeros below the critical point
        # and no zeros above the critical point. These zeros correspond to the spinodal.
        # As the temperature increases, the points on the spinodal approach the critical point, allowing
        # us to estimate it.
        if regime == 'low':
            cur_temp = parameter_list[0] - 50
            cur_crit_phi = False
            inc = 10
            while inc > .01:    
                theory_class = theory_module.Theory_Class(parameter_list, cur_temp)
                cur_roots = theory_class.GCT_ddf.roots()
                if len(cur_roots) == 2:
                    cur_crit_T = cur_temp
                    cur_crit_phi = np.average(cur_roots)
                    cur_temp += inc
                elif cur_crit_phi:
                    inc /= 10
                    cur_temp -= inc * 9
                else:
                    cur_temp -= 50
            crit_phi = cur_crit_phi
            crit_var = cur_crit_T
        else:
            cur_temp = parameter_list[0] + 50
            cur_crit_phi = False
            inc = 10
            while inc > .01:    
                theory_class = theory_module.Theory_Class(parameter_list, cur_temp)
                cur_roots = theory_class.GCT_ddf.roots()
                if len(cur_roots) == 2:
                    cur_crit_T = cur_temp
                    cur_crit_phi = np.average(cur_roots)
                    cur_temp -= inc
                elif cur_crit_phi:
                    inc /= 10
                    cur_temp += inc * 9
                else:
                    cur_temp += 50
            crit_phi = cur_crit_phi
            crit_var = cur_crit_T
    else:
        theory_class = theory_module.Theory_Class(parameter_list)
        # Function that builds the critical_phi_function. This internal function  finds the value of 
        # crit_phi by solving for when dddf is zero. This value is independent of T.
        crit_phi_fx = theory_class.build_critical_phi_function()    
        crit_phi = optimize.root(crit_phi_fx, 0.01).x[0]
        # Find the value of crit_T by explicitly solving for when ddf is zero
        # This value is dependent on phi, so we substitute in the value for crit_phi
        crit_var_fx = theory_class.build_spinodal_function_var(crit_phi)
        crit_var = optimize.root(crit_var_fx, 300).x[0]
    # Modify dilute_array and dense_array in case the value of crit_phi is larger than a value in dense_array or smaller than a value in dilute_array
    dilute_array, dense_array = array_fixer(dilute_array, dense_array, crit_phi)
    
    # Call partial_binodal_loop to do the heavy lifting and actually calculate the relevant binodal points
    fake_dilute_array = np.pad(dilute_array, [(0, 0), (0, 1)], mode='constant', constant_values=0)
    fake_dense_array = np.pad(dense_array, [(0, 0), (0, 1)], mode='constant', constant_values=1)
    data_array = np.concatenate((fake_dilute_array, fake_dense_array), axis=0)
    data_array = data_array[data_array[:,1].argsort()[::-1]]
    full_dict = {}
    if mode == 'GCT':
        full_dict[regime] = partial_binodal_loop(mode, regime, data_array, crit_phi, crit_var, temperature_offset, parameter_list = parameter_list)
    else:
        full_dict[regime] = partial_binodal_loop(mode, regime, data_array, crit_phi, crit_var, temperature_offset, theory_class = theory_class)        

    #===========================================================================
    ### Future proofing for the ability to add LCST data
    # if array_fit_dict["low"] == 1:
    #     full_dict['low_dilute'] = partial_binodal_loop(theory_class, 'low', low_dilute_array, 'dilute', crit_phi, crit_var, temperature_offset)
    #     full_dict['low_dense'] = partial_binodal_loop(theory_class, 'low', low_dense_array, 'dense', crit_phi, crit_var, temperature_offset)
    # if array_fit_dict["high"] == 1:
    #     full_dict['high_dilute'] = partial_binodal_loop(theory_class, 'high', high_dilute_array, 'dilute', crit_phi, crit_var, temperature_offset)
    #     full_dict['high_dense'] = partial_binodal_loop(theory_class, 'high', high_dense_array, 'dense', crit_phi, crit_var, temperature_offset)
    #===========================================================================
    
    # Here we calculate the total normalized residuals and print the final value
    # This allows us to see how the normalized residuals evolve over the course of the optimization
    full_length = 0
    full_residuals = 0
    full_phi_array = []
    full_var_array = []
    full_error_array = []
    for key in full_dict:
        full_length += full_dict[key]['length']
        full_residuals += full_dict[key]['residuals']
        full_phi_array.extend(full_dict[key]['phi'])
        full_var_array.extend(full_dict[key]['var'])
        full_error_array.extend(full_dict[key]['error'])
    full_residuals /= full_length
    print("Normalized Residuals: %s" %(str(10 ** np.sqrt(full_residuals))))
    if error:
        print("Error List: " + ', '.join([str(err) for err in full_error_array]))
    
    # plotter is an argument for this function
    # The partial binodal is plotted along with the data if plotter == 1
    if plotter:

        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        # local imports are not (in general) recommended, but this means that we only import these libraries IF 
        # we try and plot. This has the advnatage that we can run this code on a headless server where (lack of) GTK backends 
        # might make importing mpl or sns hard
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # The y-variable for our phase diagrams
        y_variable = config.y_variable 
        
        # The y-axis label for our diagrams
        y_label = config.y_label 
        
        # We create pandas dataframes from the partial binodal arrays and our data
        binodal_df = pd.DataFrame({'Volume Fraction': full_phi_array, y_variable: full_var_array})
        our_df = pd.DataFrame(np.append(dilute_array - [0, temperature_offset, 0], dense_array - [0, temperature_offset, 0], axis=0), columns = ['Volume Fraction', y_variable, 'Volume Fraction Error'])
        
        #=======================================================================
        ### Future proofing for UCST + LCST data
        # if array_fit_dict["low"] == 1 and array_fit_dict["high"] == 0:
        #     our_df = pd.DataFrame(np.append(low_dilute_array - [0, temperature_offset], low_dense_array - [0, temperature_offset], axis=0), columns = ['Volume_Fraction', y_variable])
        # if array_fit_dict["low"] == 0 and array_fit_dict["high"] == 1:
        #     our_df = pd.DataFrame(np.append(high_dilute_array - [0, temperature_offset], high_dense_array - [0, temperature_offset], axis=0), columns = ['Volume_Fraction', y_variable])
        # if array_fit_dict["low"] == 1 and array_fit_dict["high"] == 1:
        #     our_df = pd.DataFrame(np.append(low_dilute_array - [0, temperature_offset], low_dense_array - [0, temperature_offset], high_dilute_array - [0, temperature_offset], high_dense_array - [0, temperature_offset], axis=0), columns = ['Volume_Fraction', y_variable])
        #=======================================================================
        
        # We graph and show the data with the binodal
        fig = plt.figure(figsize=(12,5))
        gs = fig.add_gridspec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        
        fig.suptitle(title)

        ax1.set(xlabel='Volume Fraction', ylabel = y_label)
        ax2.set(xlabel='Volume Fraction', ylabel = y_label)
        
        ax1.errorbar(x = binodal_df['Volume Fraction'], y = binodal_df[y_variable], fmt='.', color = sns.color_palette()[0], capsize=2, markersize=10, markeredgewidth=0.5, markeredgecolor='white', label='binodal')
        if our_df['Volume Fraction Error'].isnull().all():
            ax1.errorbar(x = our_df['Volume Fraction'], y = our_df[y_variable], fmt='^', color=sns.color_palette()[1], ecolor='black', capsize=2, markersize=8, markeredgewidth=0.5, markeredgecolor='white', label='data')
        else:
            ax1.errorbar(x = our_df['Volume Fraction'], y = our_df[y_variable], xerr = our_df['Volume Fraction Error'], fmt='^', color=sns.color_palette()[1], ecolor='black', capsize=2, markersize=8, markeredgewidth=0.5, markeredgecolor='white', label='data')
        ax1.legend(loc='upper right')
        
        ax2.set(xscale="log")
        ax2.set(xlim=(.000001, 1))
        ax2.errorbar(x = binodal_df['Volume Fraction'], y = binodal_df[y_variable], fmt='.', color = sns.color_palette()[0], capsize=2, markersize=10, markeredgewidth=0.5, markeredgecolor='white', label='binodal')
        if our_df['Volume Fraction Error'].isnull().all():
            ax2.errorbar(x = our_df['Volume Fraction'], y = our_df[y_variable], fmt='^', color=sns.color_palette()[1], ecolor='black', capsize=2, markersize=8, markeredgewidth=0.5, markeredgecolor='white', label='data')
        else:
            ax2.errorbar(x = our_df['Volume Fraction'], y = our_df[y_variable], xerr = our_df['Volume Fraction Error'], fmt='^', color=sns.color_palette()[1], ecolor='black', capsize=2, markersize=8, markeredgewidth=0.5, markeredgecolor='white', label='data')
        ax2.legend(loc='upper left')
        
        plt.tight_layout()
        plt.show()
        
    # The function returns the normalized residuals regardless of whether the data were plotted
    return full_residuals


