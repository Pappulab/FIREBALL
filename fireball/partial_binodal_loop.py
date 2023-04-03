'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

# import required packages
import importlib
import numpy as np
from scipy import optimize

from .fireballexceptions import FireballException
from . import config


def partial_binodal_loop(mode, regime, data_array, crit_phi, crit_var,
                         temperature_offset, parameter_list = None,
                         theory_class = None):
    """

    This function is called by partial_binodal to loop through dilute and dense arm
    measurements and calculate the binodal at the given concentrations, given the passed
    parameters. The binodal values are used to calculate residuals from the data.
    The data are all returned in a dictionary. As opposed to partial_binodal_loop,
    this function calculates the binodal as is done in full_binodal

    Parameters
    ------------
    mode : str
        The name of the theory module to use.

    regime : str (either "low" or "high")
        If "low", then the algorithm assumes an upper critical system
        If "high", then the algorithm assumes a lower critical system
        By default, fireball-fit and fireball-draw use "low"

    data_array : np.ndarry (3xm)
        A 4D array where column 0 is the volume fraction, column 1 is the temperature,
        column 2 is the volume fraction error, and column 3 is the arm type (dilute=0, dense=1).

    crit_phi : float
        The critical volume fraction. Should be determined directly in partial_binodal.py
    
    crit_var : float
        The critical value for the independent variable (salt, temp, etc).
        Should be determined directly in partial_binodal.py

    temperature_offset : float
        Converts the y-axis units when plotting the data. Should be set to 273.15 to convert the units
        from Kelvin to degrees Celsius. Should bet set to 0 to keep the units in Kelvin.

    parameter_list : list
        List that defines the set of parameters used to generate the partial binodal for the given theory module.
        Required for GCT calculations. theory_class should instead be used for FH calculations.
        
    theory_class : class
        The theory module class initialized in the partial_binodal.py module.
        Required for FH calculations. parameter_list should instead be used for GCT calculations.

    Returns
    --------
    dictionary
       A dictionary that describes the data in the computed binodal.
       The dictionary has 5 keys - phi, var, residuals, error, and length, which correspond to
       the volume concentrations (np array), the values of the independent variable (temp, salt, etc; np array), the 
       rescaled residuals of the fit from the data (float), the computational error in computing the binodal 
       points that correspond to the data (np array), and the number of data points (float) that were used to calculate 
       residuals, respectively

    """
    
    # Load in the required information depending on which theory is being used
    if mode == 'GCT':
        if parameter_list == None:
            raise FireballException('parameter_list in partial_binodal_loop must not be False for GCT calculations')
        theory_module = importlib.import_module('fireball.theory.' + mode)
    else:
        if theory_class == None:
            raise FireballException('theory_class in partial_binodal_loop must not be False for FH calculations')
    
    if regime == 'low':
        regime = 0
    elif regime == 'high':
        regime = 1

    # Initialize the arrays we will need as we go through the data
    phi_array = np.array([])
    var_array = np.array([])
    error_set = []
    residuals = 0
    relative_size = 0
    
    # Initialize the remaining parameters we will need as we construct the partial binodal
    # increment determines the gap between consecutive points on the binodal. Smaller values
    # of increment mean more points will be calculated while trying to calculate the residuals.
    # This helps the binodal to be more accurate, but will slow down the calculations. We want
    # to choose the largest possible value of increment, while maintaining accuracy.
    i = 0
    current_var = crit_var
    increment = config.partial_binodal_increment
    dilute_high_bound = 1
    dense_low_bound = 0
    binodal_check = 0
    
    # This loop generates partial binodal and calculates the residuals along the way
    while i < (data_array.shape[0]):
    
    # In case the independent variable values from the data are very far from crit_var
    # or from each other, this loop constructs intermediate points along the binodal.
    # This gives us stricter bounds for the next phi values from the data.
    # The smaller increment is, the more intermediate points will be calculated.
    # This makes the later calculations more accurate, but takes longer.
    # We want to use the largest value for increment we can that still gives accurate calculations.
    # Residuals are only calculated when the calculated binodal points correspond to
    # independent variable values from the data, not for the intermediate points 
        if regime == 0:
            # UCST system
            
            if data_array[i][1] > current_var:
                # The data points are above the UCST, so we calculate their distance from
                # the critical point to measure pseudo-residuals.
                residuals += (np.log10(data_array[i][0] / crit_phi)) ** 2
                residuals += ((data_array[i][1] - crit_var) / 10) ** 2
                relative_size += 1
                i += 1
                continue
            
            if data_array[i][1] < (current_var - increment):
                # The data points are too far from the last set of calculated values,
                # so we calculate intermediate values.
                binodal_check = 0
                current_var = current_var - increment
                
            elif binodal_check == 1 and current_var == data_array[i][1]:
                # We already calculated this binodal point in the last loop but for the other
                # binodal arm! So we can directly calculate residuals based on those values.
                if data_array[i][3] == 0:
                    residuals += (np.log10(data_array[i][0] / phi_L)) ** 2
                elif data_array[i][3] == 1:
                    residuals += (np.log10(data_array[i][0] / phi_H)) ** 2
                    #residuals += (data_array[i][0] - phi_H) ** 2
                else:
                    raise FireballException('Partial binodal can not read correct arm types')
                relative_size += 1
                error_set.append(binodal.fun)
                i += 1
                continue
            
            else:
                # We are within range to calculate binodal points that correspond to our real
                # data points, so we will go ahead and do that
                binodal_check = 1
                current_var = data_array[i][1]
        else:
            # LCST system
            
            if data_array[i][1] < current_var:
                # The data points are below the LCST, so we calculate their distance from
                # the critical point to measure pseudo-residuals.
                residuals += (np.log10(data_array[i][0] / crit_phi)) ** 2
                residuals += ((data_array[i][1] - crit_var) / 10) ** 2
                relative_size += 1
                i += 1
                continue
            
            if data_array[i][1] > (current_var - increment):
                # The data points are too far from the last set of calculated values,
                # so we calculate intermediate values.
                binodal_check = 0
                current_var = current_var + increment
                
            elif binodal_check == 1 and current_var == data_array[i][1]:
                # We already calculated this binodal point in the last loop but for the other
                # binodal arm! So we can directly calculate residuals based on those values.
                if data_array[i][3] == 0:
                    residuals += (np.log10(data_array[i][0] / phi_L)) ** 2
                elif data_array[i][3] == 1:
                    residuals += (np.log10(data_array[i][0] / phi_H)) ** 2
                    #residuals += (data_array[i][0] - phi_H) ** 2
                else:
                    raise FireballException('Partial binodal can not read correct arm types')
                relative_size += 1
                error_set.append(binodal.fun)
                i += 1
                continue
                      
            else:
                # We are within range to calculate binodal points that correspond to our real
                # data points, so we will go ahead and do that
                binodal_check = 1
                current_var = data_array[i][1]

        if mode == 'GCT':
            # Instantiate a theory_class based on the parameter_list and the current temperature
            # Put together the values and functions needed to calculate the binodal point
            theory_class = theory_module.Theory_Class(parameter_list, current_var)
            spinodal_function_phi = theory_class.build_spinodal_function_phi()
            binodal_function = theory_class.build_binodal_function_no_phi()
            jacobian_function = theory_class.build_jacobian_function_no_phi()
            spinodal_low_phi = optimize.brentq(spinodal_function_phi, 10 ** -30, crit_phi - 0.0000001)
            spinodal_high_phi = optimize.brentq(spinodal_function_phi, crit_phi + 0.0000001, theory_class.cur_max)
            dilute_high_bound = np.min([dilute_high_bound, spinodal_low_phi])
            dense_low_bound = np.max([dense_low_bound, spinodal_high_phi])
            dense_high_bound_1 = (-1 * theory_class.GCT_T(current_var) / parameter_list[2] - 16 / 3) / (theory_class.F * np.sqrt(theory_class.N))
            dense_high_bound_2 = dense_low_bound * 2
            dense_high_bound_3 = theory_class.cur_max
            dense_high_bound = np.min([dense_high_bound_1, dense_high_bound_2, dense_high_bound_3])
            bounds = ((dense_low_bound, dense_high_bound), (10 ** -12, dilute_high_bound))
            x0 = [(dense_low_bound + 1) / 5, dilute_high_bound * 0.9]
        else:
            # Put together the values and functions needed to calculate the binodal point
            spinodal_function_phi = theory_class.build_spinodal_function_phi(current_var)
            binodal_function = theory_class.build_binodal_function_no_phi(current_var)
            jacobian_function = theory_class.build_jacobian_function_no_phi(current_var)
            spinodal_low_phi = optimize.brentq(spinodal_function_phi, 10 ** -30, crit_phi - 0.0000001)
            spinodal_high_phi = optimize.brentq(spinodal_function_phi, crit_phi + 0.0000001, 0.99999)
            dilute_high_bound = np.min([dilute_high_bound, spinodal_low_phi])
            dense_low_bound = np.max([dense_low_bound, spinodal_high_phi])
            bounds = ((dense_low_bound, .9999999), (10 ** -10, dilute_high_bound))
            x0 = [(dense_low_bound + 1) / 5, dilute_high_bound * 0.9]

        tol = 100 ** -10
        method = "L-BFGS-B"

        # Perform the minimization
        binodal = optimize.minimize(binodal_function, x0, bounds = bounds, method = method, jac = jacobian_function, tol = tol)
        
        # Define new bounds for the next iteration based on the results of our last minimization
        dense_low_bound = phi_H = binodal.x[0]
        dilute_high_bound = phi_L = binodal.x[1]
                 
        # Append binodal points to the arrays
        phi_array = np.append(phi_array, [phi_L, phi_H])
        var_array = np.append(var_array, [current_var - temperature_offset, current_var - temperature_offset])
        #error_set.append(binodal.fun)
        if binodal_check == 1:
            if data_array[i][3] == 0:
                residuals += (np.log10(data_array[i][0] / phi_L)) ** 2
            elif data_array[i][3] == 1:
                residuals += (np.log10(data_array[i][0] / phi_H)) ** 2
                #residuals += (data_array[i][0] - phi_H) ** 2
            else:
                print(data_array)
                raise FireballException('Partial binodal can not read correct arm types')
            relative_size += 1
            error_set.append(binodal.fun)
            i += 1
    
    # Return the final dictionary with all required values
    array_dict = {'phi': phi_array, 'var': var_array, 'residuals': residuals, 'error': error_set, 'length': relative_size}
    return array_dict


