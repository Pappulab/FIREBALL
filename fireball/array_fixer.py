'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

# import required packages
import numpy as np

from fireball.fireballexceptions import FireballException

# The next few lines modify dilute_array and dense_array in case the value of crit_phi is
# larger than a value in dense_array or smaller than a value in dilute_array
def array_fixer(dilute_array, dense_array, crit_phi):

    """
    Function that modifies the dilute and dense arrays in case the value of crit_phi is larger than
    a value in dense_array or smaller than a value in dilute_array.    

    Parameters
    -------------
    dilute_array : np.ndarray (2xm)
        A 2D array where column 0 is the dilute phase concentration and column 1 is 
        the independent variable (temp, salt, etc) of the binodal.

    dense_array : np.ndarray (2xm)
        A 2D array where column 0 is the dense phase concentration and column 1 is 
        the independent variable (temp, salt, etc) of the binodal.
        
    crit_phi : float
        The critical phi value determined from setting dddf = 0 and solving for phi.

    Returns
    ----------

    tuple of 2 np.ndarrays (2xm)
        Returns the new dilute_array and dense_array

    """

    if dilute_array.shape[0] != 0:
        while crit_phi < dilute_array[0][0]:
            dense_array = np.insert(dense_array, 0, dilute_array[0], axis = 0)
            dilute_array = np.delete(dilute_array, 0, axis = 0)
            if dilute_array.shape[0] == 0:
                break
            
    if dense_array.shape[0] != 0:
        while crit_phi > dense_array[0][0]:
            dilute_array = np.insert(dilute_array, 0, dense_array[0], axis = 0)
            dense_array = np.delete(dense_array, 0, axis = 0)
            if dense_array.shape[0] == 0:
                break
    
    return (dilute_array, dense_array)


