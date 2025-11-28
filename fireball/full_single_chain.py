'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

# Import required packages
import importlib
import numpy as np
import pandas as pd
import math
from scipy import optimize

from fireball.fireballexceptions import FireballException
from fireball import config

def build_theoretical_coil_globule(parameter_list, mode):
    """
    Simple function that returns a dataframe with temperature and Rg. If the goal is simply to examine
    how coil-globule transitions change as a function of the input parameters, this is the easiest way to do this, although 
    fundamentally this function is just a wrapper around fireball.single_chain.coil_globule_maker().

    Parameters
    -----------    
    parameter_list : list
        list of parameters as defined in the given theory module (default is GCT_Single_Chain)

    mode : str
        The name of the theory module to use. By default, fireball-single-chain uses GCT_Single_Chain
        
    Returns
    ---------
    (pd.dataframe, list, list)

        Returns a tuple of three elements:

        element 0: dataframe with the theoretical transition. 
        element 1: an ordered list of temperatures that map to the coil-globule transition
        element 2: an ordered list of reduced temperatures that map to the coil-globule transition
        element 3: an ordered list of Rg values that map to the coil-globule transition
        element 4: an ordered list of alpha values that map to the coil-globule transition

    """

    return coil_globule_maker(parameter_list, mode, None, outname=None, plotter=False)


def coil_globule_maker(parameter_list, mode, data_array, temperature_offset=None, increment=None, outname=None, plotter=True):
    """
    Function for generating coil_globule transitions using the given theory method.

    Parameters
    --------------
    parameter_list : list
        list of parameters as defined in the given theory module (default is GCT_Single_Chain)

    mode : str
        The name of the theory module to use. By default, fireball-single-chain-fit and fireball-single-chain-draw use GCT_Single_Chain
        
    data_array : np.ndarray (2D)
        Array of points along the coil-globule transition. Includes both Rg and temperature.

    temperature_offset : float
        Converts the y-axis units when plotting the data. Should be set to 273.15 to convert the units
        from Kelvin to degrees Celsius. Should bet set to 0 to keep the units in Kelvin.

    increment : float
        The increment between consecutive values of the independent variable (temp, salt, etc.)

    outname : str or None
        Name of output file that is written containing binodal data. Note if set to None no file is written. Default = None

    plotter : bool
        Flag which, if set to true, ensures that the resulting binodal are plotted to screen and saved to a local file (binodal.pdf
        and binodal_log.pdf). Default = True

    Returns
    ---------
    (pd.dataframe, list, list, list)

        Returns a tuple of four elements:

        element 0: dataframe with the theoretical transition. 
        element 1: an ordered list of temperatures that map to the coil-globule transition
        element 2: an ordered list of reduced temperatures that map to the coil-globule transition
        element 3: an ordered list of Rg values that map to the coil-globule transition
        element 4: an ordered list of alpha values that map to the coil-globule transition

    """
    
    # The theory we are currently using. We will instantiate a Class based on this theory.
    theory_module = importlib.import_module('fireball.theory.' + mode)
    theory_class = theory_module.Theory_Class(parameter_list)

    # temperature_offset converts the y-axis units when plotting the data.
    # temperature_offset should bet set to 273.15 to convert the units from Kelvin to degrees Celsius
    # temperature_offset should bet set to 0 to keep the units in Kelvin
    if temperature_offset is None:
        temperature_offset = config.temperature_offset
    x_label = config.y_label
    
    # Directly calculate all of the temperatures and Rg values
    real_temp_array = [t_val - temperature_offset for t_val in theory_class.real_t_vals]
    reduced_temp_array = theory_class.reduced_t_vals
    rg_squared_array = theory_class.rg_vals
    rg_array = [np.sqrt(rg2) for rg2 in theory_class.rg_vals]
    alpha_array = theory_class.alpha_vals
    
    # Create pandas dataframes from the coil-globule transition arrays and our data
    try:
        coil_globule_df = pd.DataFrame({'Temperature': real_temp_array, 'Reduced Temperature': reduced_temp_array, 'Rg': rg_array, 'alpha': alpha_array})
    except ValueError:
        coil_globule_df = pd.DataFrame({'Temperature': [], 'Reduced Temperature': [], 'Rg': [], 'alpha': []})

    if outname is not None:
        # Write out the coil-globule transition to a file
        with open(outname,'w') as fh:
            for i in range(len(real_temp_array)):
                fh.write('%.3f, %.3f, %.3f, %.3f\n' % (real_temp_array[i], reduced_temp_array[i], rg_array[i], alpha_array[i]))

    if data_array is not None:
        our_df = pd.DataFrame(data_array - [temperature_offset, 0, 0], columns = ['Temperature', 'Rg', 'Rg Error'])
        # Concatenate the dataframes and create a new column describing if the data is calculated or measured
        full_frame = pd.concat([coil_globule_df.assign(dataset = 'Fit'), our_df.assign(dataset = 'data')])        
    else:
        full_frame = coil_globule_df.assign(dataset = 'Fit')

    # if no plotting just return
    if plotter is False:
        return (coil_globule_df, real_temp_array, reduced_temp_array, rg_array, alpha_array)

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    # local imports are not (in general) recommended, but this means that we only import these libraries IF 
    # we try and plot. This has the advantage that we can run this code on a headless server where (lack of) GTK backends 
    # might make importing mpl or sns hard
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Graph, save, and show the data with the calculated coil-globule transition
    # ax1 shows the data on a linear scale
    # ax2 shows the data on a semi-log scale
    fig = plt.figure(figsize=(6, 6))
    gs = fig.add_gridspec(1, 1)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set(xlabel=x_label, ylabel = 'Rg (nm)')

    ax1.errorbar(x = coil_globule_df['Temperature'], y = coil_globule_df['Rg'], fmt='-', color = sns.color_palette()[0], capsize=2, markersize=10, markeredgewidth=0.5, markeredgecolor='white', label='Fit')
    if data_array is not None:
        if our_df['Rg Error'].isnull().all():
            ax1.errorbar(x = our_df['Temperature'], y = our_df['Rg'], fmt='^', color=sns.color_palette()[1], ecolor='black', capsize=2, markersize=8, markeredgewidth=0.5, markeredgecolor='white', label='Data')
        else:
            ax1.errorbar(x = our_df['Temperature'], y = our_df['Rg'], xerr = our_df['Rg Error'], fmt='^', color=sns.color_palette()[1], ecolor='black', capsize=2, markersize=8, markeredgewidth=0.5, markeredgecolor='white', label='Data')        
    ax1.legend(loc='upper left')
    
    #===========================================================================
    # ## For comparisons with Raos and Allegra figures
    # ax2 = fig.add_subplot(gs[0, 1])
    # ax2.set(xlabel='Reduced temperature', ylabel = 'alpha')
    # ax2.errorbar(x = coil_globule_df['Reduced Temperature'], y = coil_globule_df['alpha'], fmt='-', color = sns.color_palette()[0], capsize=2, markersize=10, markeredgewidth=0.5, markeredgecolor='white', label='Fit')
    # ax2.set_xlim(-0.7, 0)
    # ax2.set_ylim(0, 1)
    # ax2.legend(loc='upper right')
    #===========================================================================
    
    plt.tight_layout()
    plt.savefig('coil_globule.pdf', bbox_inches='tight')
    plt.show()

    return (coil_globule_df, real_temp_array, reduced_temp_array, rg_array, alpha_array)
    
