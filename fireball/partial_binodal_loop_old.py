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


def partial_binodal_loop(theory_class, regime, arm_array, arm_type, crit_phi, crit_var, temperature_offset):
    """

    This function is called by partial_binodal to loop through dilute and dense arm measurements and calculate the binodal at the given 
    concentrations, given the passed parameters. The binodal values are used to calculate residuals from the data. The data are all returned 
    in a dictionary.

    Parameters
    ------------
    theory_class : class
        The theory module class initialized in the partial_binodal.py module. By default, fireball-fit and fireball-draw use flory_huggins_3B
        
    regime : str (either "low" or "high")
        If "low", then the algorithm assumes an upper critical system
        If "high", then the algorithm assumes a lower critical system
        By default, fireball-fit and fireball-draw use "low"

    arm_array : np.ndarry (2xm)
        A 2D array where column 0 is the volume fraction and column 1 is the temperature.
        This array should belong to either the dilute or dense arm
    
    arm_type : str (either "dense" or "dilute")
        Describes whether the data in arm_array belong to a dilute or dense arm. This tells the algorithm how to sift through the data

    crit_phi : float
        The critical volume fraction. Should be determined directly in partial_binodal.py
    
    crit_var : float
        The critical value for the independent variable (salt, temp, etc).
        Should be determined directly in partial_binodal.py

    temperature_offset : float
        Value which is used as a global offset to define the region plotted - i.e. this is the value that corresponds to 0 on the 
        temperature scale being used. Recommend 273.15 if data are reported in Kelvin.

    Returns
    --------
    dictionary
       A dictionary that describes the data in the computed binodal.
       The dictionary has 5 keys - phi, var, residuals, error, and length, which correspond to
       the volume concentrations (np array), the values of the independent variable (temp, salt, etc; np array), the 
       sum of the square residuals of the fit from the dat (float), the computational error in computing the binodal 
       points that correspond to the data (np array), and the number of data points (float) that were used to calculate 
       residuals, respectively

    """
    
    if arm_type == 'dilute':
        arm_type = 0
    elif arm_type == 'dense':
        arm_type = 1
    if regime == 'low':
        regime = 0
    elif regime == 'high':
        regime = 1

    # We initialize the arrays we will need as we go through the data
    phi_array = np.array([])
    var_array = np.array([])
    error_set = []
    residuals = 0
    relative_size = 0
    
    # low_crit_var and high_crit_var are temporary placeholders for future-proofing for UCST + LCST data
    low_crit_var = crit_var
    high_crit_var = crit_var

    # Here we initialize the remaining parameters we will need as we construct the partial binodal 
    i = 0
    current_phi = crit_phi
    var_high_bound = low_crit_var
    var_low_bound = high_crit_var
    
    # As we construct the low partial binodal, we begin close to crit_phi and move outwards
    # Here we determine a bound for only the first phi value along each arm. The bounds for the following phi values are determined by the prior phi values 
    # These values may need to be changed depending on the system
    # Bounds that are further from the critical point are more reliable, but can only be used if the data do not have points near the critical point
    # We want these points to be as large as possible while still being valid for the first point along each arm, with the understanding that the dilute data determines the dense arm and vice-versa        
    dilute_high_bound = crit_phi - config.dilute_first_phi
    dense_low_bound = crit_phi + config.dense_first_phi

    # This loop generates the opposite arm of the partial binodal and calculates the residuals along the way
    while i < (arm_array.shape[0]):
    
    # In case the phi values from the data are very far from crit_phi or from each other, this loop constructs intermediate points along the dense arm
    # This gives us stricter bounds for the next phi value from the data
    # The smaller dilute_intermediate_increment is, the more intermediate points will be calculated
    # This makes the later calculations more accurate, but takes longer
    # We want to use the largest value for dilute_intermediate_increment we can that still gives accurate calculations
    # current_phi is the phi value along the dilute arm
    # If the next phi value from the data is close enough to the last binodal point found, then we perform a minimization using the value from the data
    # In this case we included the residual and the error at the end
    
        if arm_type == 0:
            if arm_array[i][0] > current_phi:
                residuals += (arm_array[i][1] - low_crit_var) ** 2
                relative_size += 1
                i += 1
                continue
            if arm_array[i][0] < (current_phi - current_phi * config.dilute_intermediate_increment):
                binodal_check = 0
                current_phi = current_phi - current_phi * config.dilute_intermediate_increment
            else:
                binodal_check = 1
                current_phi = arm_array[i][0]
        elif arm_type == 1:
            if arm_array[i][0] < current_phi:
                residuals += (arm_array[i][1] - low_crit_var) ** 2
                relative_size += 1
                i += 1
                continue
            if arm_array[i][0] > (current_phi + config.dense_intermediate_increment):
                binodal_check = 0
                current_phi = current_phi + config.dense_intermediate_increment
            else:
                binodal_check = 1
                current_phi = arm_array[i][0]                
                
        # var_low_bound/var_high_bound is the temperature at the corresponding point on the spinodal
        spinodal_function_var = theory_class.build_spinodal_function_var(current_phi)
        
        if regime == 0:
            #var_low_bound = max(var_high_bound - 200, 0.000001)
            #var_low_bound = max(var_high_bound - 2, optimize.brentq(spinodal_function_var, .00001, var_high_bound))
            var_low_bound = optimize.brentq(spinodal_function_var, .00001, var_high_bound)
            
        elif regime == 1:
            #var_high_bound = var_low_bound + 100
            var_high_bound = optimize.brentq(spinodal_function_var, var_low_bound, var_low_bound * 10000)
            
        # functions binodal_function and jacobian are used to find the desired point on the binodal 
        binodal_function = theory_class.build_binodal_function(current_phi)
        jacobian_function = theory_class.build_jacobian_function(current_phi)

        # Here we define our minimization parameters to find the desired point on the binodal            
        if arm_type == 0 and regime == 0:
            bounds = ((dense_low_bound, .99), (var_low_bound, var_high_bound))
            x0 = [(dense_low_bound * 2), var_high_bound - 1]
        elif arm_type == 0 and regime == 1:
            bounds = ((dense_low_bound, .99), (var_low_bound, var_high_bound))
            x0 = [(dense_low_bound * 2), (var_low_bound + 1)]
        elif arm_type == 1 and regime == 0:
            bounds = ((0.000000000001, dilute_high_bound), (var_low_bound, var_high_bound))
            x0 = [dilute_high_bound * 0.8, var_high_bound - 1]
        elif arm_type == 1 and regime == 1:
            bounds = ((0.000000000001, dilute_high_bound), (var_low_bound, var_high_bound))
            x0 = [dilute_high_bound * 0.95, var_low_bound + 1]
        tol = 100 ** -10
        method = "L-BFGS-B"
        
        # Here we perform the minimization
        binodal = optimize.minimize(binodal_function, x0, bounds = bounds, method = method, jac = jacobian_function, tol = tol)
        
        # Here we define new bounds for the next iteration based on the results of our last minimization
        if arm_type == 0:
            dense_low_bound = new_phi = binodal.x[0]
        elif arm_type == 1:
            dilute_high_bound = new_phi = binodal.x[0]
        if regime == 0:
            var_high_bound = new_var = binodal.x[1]
        elif regime == 1:
            var_low_bound = new_var = binodal.x[1]
         
        # Here we append our binodal points to our arrays
        phi_array = np.append(phi_array, [current_phi, new_phi])
        var_array = np.append(var_array, [new_var - temperature_offset, new_var - temperature_offset])
        #error_set.append(binodal.fun)
        if binodal_check == 1:
            residuals += (arm_array[i][1] - new_var) ** 2
            relative_size += 1
            error_set.append(binodal.fun)
            i += 1
    
    array_dict = {'phi': phi_array, 'var': var_array, 'residuals': residuals, 'error': error_set, 'length': relative_size}
    return array_dict


