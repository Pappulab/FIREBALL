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
from . import config

def residual_finder(free_parameter_list, fixed_parameter_list, free_parameter_checklist,
                    mode, data_array, temperature_offset, plotter=False, error=False, title=''):
    """

    This program is used to construct a partial coil-globule transition using a given set
    of parameters and a theory module.
    
    This program finds the points along the coil-globule transition that correspond to
    the temperature values of the data. The residuals for these data points are the difference
    in Rg values between the data points and the points on the partial coil-globule
    transition. This program can be used in conjunction with optimizer_single_chain.py
    in order to find the best fit parameters.

    This function is written such that it is compatible with scipy.optimize.minimize,
    i.e. first argument is a vector of parameters to be optimized over. This function is not
    really much use OTHER than in the context of being passed to minimize, which happens 
    in the optimizer_single_chain.py module.

    Parameters
    ------------
    free_parameter_list : list
        List that defines the set of parameters set to be optimized for the given theory
        module.

    fixed_parameter_list : list
        List that defines the set of parameters set to be fixed for the given theory module.

    free_parameter_checklist : list
        List that defines which parameters should be fixed and which should be optimized.
    
    mode : str
        The name of the theory module to use.
        
    data_array : np.ndarray (3xm)
        A 2D array where column 0 is the temperature, column 1 is the Rg, and column 3 is
        the Rg error. This is automatically generated as output element 0 from the
        fireball.io_single_chain.read_file() function

    temperature_offset : float
        Converts the y-axis units when plotting the data. Should be set to 273.15 to convert the units
        from Kelvin to degrees Celsius. Should bet set to 0 to keep the units in Kelvin.

    plotter : bool
        Flag which, if set to true, means that the partial coil-globule transition fit
        is plotted to screen (only). Useful when debugging a fit, but in general
        not necessary. Default = False.

    error : bool
        Flag which, if set to true, prints out the list of errors in calculating the
        coil-globule to globule transition points that correspond to given data points.
        This can be helpful in troubleshooting whether the algorithm is correctly
        computing points.
    
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
    
    # Create an instantiation of theory class with the desired parameters and initialize
    # the required variables
    theory_class = theory_module.Theory_Class(parameter_list, partial_check=True)
    temp_array = []
    rg_array = []
    full_length = 0
    full_residuals = 0
    
    # Here we calculate the sum of the square-residuals and print the final value
    # This allows us to see how the sum of the square-residuals evolves over the course of the optimization
    for data_point in data_array:
        if data_point[0] < theory_class.GCT_T_inverse(theory_class.cur_min):
            # If the temperature is below the extrapolation threshold for the compression
            # factor, we calculate the Rg based on the minimum allowable temperature,
            # given the parameters used.
            rg = np.sqrt(theory_class.GCT_rg2(theory_class.q_interp(theory_class.cur_min)))
        elif theory_class.GCT_T_inverse(theory_class.cur_max) < data_point[0]:            
            # If the temperature is above the extrapolation threshold for the compression
            # factor, we calculate the Rg based on a linear extrapolation of the Rg plot.
            temp_2 = theory_class.cur_max
            temp_1 = theory_class.GCT_T(theory_class.GCT_T_inverse(temp_2) - 1)
            rg_1 = np.sqrt(theory_class.GCT_rg2(theory_class.q_interp(temp_1)))
            rg_2 = np.sqrt(theory_class.GCT_rg2(theory_class.q_interp(temp_2)))
            rg_slope = rg_2 - rg_1
            rg = rg_2 + rg_slope * (theory_class.GCT_T_inverse(theory_class.cur_max) - data_point[0])
        else:
            # If the temperature is within the extrapolation threshold for the compression
            # factor, we directly calculate the Rg.
            rg = np.sqrt(theory_class.GCT_rg2(theory_class.q_interp(theory_class.GCT_T(data_point[0]))))
        temp_array.append(data_point[0])
        rg_array.append(rg)
        full_length += 1
        full_residuals += (rg - data_point[1]) ** 2
    full_residuals /= full_length
    print("Mean-square residuals: %s" %(str(full_residuals)))
    
    # plotter is an argument for this function
    # The partial coil-globule transition is plotted along with the data if plotter == 1
    if plotter:

        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        # local imports are not (in general) recommended, but this means that we only import these libraries IF 
        # we try and plot. This has the advnatage that we can run this code on a headless server where (lack of) GTK backends 
        # might make importing mpl or sns hard
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # The y-variable for our phase diagrams
        y_variable = 'Rg' 
        
        # The y-axis label for our diagrams
        y_label = 'Rg' 
        
        # We create pandas dataframes from the partial coil-globule transition arrays and our data
        fit_df = pd.DataFrame({'Temp': temp_array, y_variable: rg_array})
        our_df = pd.DataFrame(np.append(data_array - [temperature_offset, 0, 0], axis=0), columns = ['Temp', y_variable, 'Rg Error'])
        
        # We graph and show the data with the calculated coil-globule transition
        fig = plt.figure(figsize=(6,5))
        gs = fig.add_gridspec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        
        fig.suptitle(title)

        ax1.set(xlabel='Temperature', ylabel = y_label)
        
        ax1.errorbar(x = fit_df['Temp'], y = fit_df[y_variable], fmt='.', color = sns.color_palette()[0], capsize=2, markersize=10, markeredgewidth=0.5, markeredgecolor='white', label='Fit')
        if our_df['Rg Error'].isnull().all():
            ax1.errorbar(x = our_df['Temp'], y = our_df[y_variable], fmt='^', color=sns.color_palette()[1], ecolor='black', capsize=2, markersize=8, markeredgewidth=0.5, markeredgecolor='white', label='Data')
        else:
            ax1.errorbar(x = our_df['Temp'], y = our_df[y_variable], yerr = our_df['Rg Error'], fmt='^', color=sns.color_palette()[1], ecolor='black', capsize=2, markersize=8, markeredgewidth=0.5, markeredgecolor='white', label='Data')
        ax1.legend(loc='upper right')
                
        plt.tight_layout()
        plt.show()
        
    # The function returns the sum of the square-residuals regardless of whether the data were plotted
    return full_residuals


