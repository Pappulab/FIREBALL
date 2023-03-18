'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

import pandas as pd
import numpy as np
from . import fireballexceptions

def __validate_array(r):
    """
    Internal function that checks if the read in binodals work OK. Basically checks that binodal concentrations
    are not > 1 or < 0.

    Parameters
    ------------
    r : np.ndarray (2D)
        2D array where 0th index = concentration (and 1st = temperature)

    Returns
    ----------
    None
        No return but raises an exception if the 

    Raises
    -----------
        FireballException

    """

    if np.max(r.transpose()[0]) > 1:
        raise fireballexceptions.FireballException('Binodal concentration above 1.0. Values = %s' % (str(r)))

    if np.min(r.transpose()[0]) < 0:
        raise fireballexceptions.FireballException('Binodal concentration below 0.0. Values = %s' % (str(r)))



def read_file(filename, conv):
    """
    Function that reads in binodal data from a .xls or .csv file
        
    Parameters
    --------------
    filename : str
        Name of file to pass in. Should be a .xls or .csv file.
        Should have a 'Temp' column, a 'Dense' column, a 'Dilute' column, and 'Dense_Error' and 'Dilute_Error' columns.
        Error values can be empty if the error is unknown

    conv : float
        Conversion factor that use used in conc/CONV to convert input concentration
        into volume fraction. For mg/ml recommended factor is 1310 that converts mg/ml
        concentration into volume fraction (i.e. most tightly packed protein is 1310 mg/ml).

    Returns
    ----------
    tuple (2-place)
        Returns a 2-position tuple where 

        position = 0 is np.ndarray (2 x m) where m is number of points, col 0 is concentration
        and col 1 is temperature for the low arm of the binodal.

        position = 1 is np.ndarray (2 x n) where n is number of points, col 0 is concentration
        and col 1 is temperature for the high arm of the binodal.
        
    """

    if filename[len(filename)-3:] == 'xls':
        our_df = pd.read_excel(filename)

    elif filename[len(filename)-3:] == 'csv':
        our_df = pd.read_csv(filename, sep=r'\s*,\s*', encoding='ascii', engine='python')
    else:
        raise fireballexceptions.FireballException('\n\nERROR: Unable to parse file extension, should be .xls (excel) or .csv')
    
    our_df = our_df / [1, conv, conv, conv, conv]
    try:
        dense_array = our_df[['Dense','Temp','Dense_Error']].to_numpy()
        dilute_array = our_df[['Dilute','Temp','Dilute_Error']].to_numpy()
    except Exception as e:
        print('Unable to parse input file')
        print(e)
        raise fireballexceptions.FireballException('Error parsing file')
    
    dilute_array = dilute_array[~np.isnan(dilute_array).transpose()[0]]
    dense_array = dense_array[~np.isnan(dense_array).transpose()[0]]
    dilute_array = dilute_array[dilute_array[:,0].argsort()[::-1]]
    dense_array = dense_array[dense_array[:,0].argsort()]

    # validate arrays are OK...
    __validate_array(dilute_array)
    __validate_array(dense_array)
    
    return (dilute_array, dense_array)
