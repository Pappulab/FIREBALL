'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

import pandas as pd
import numpy as np
from . import fireballexceptions

def read_file(filename):
    """
    Function that reads in coil-globule transition data from a .xls or .csv file
    
    Parameters
    --------------
    filename : str
        Name of file to pass in. Should be a .xls or .csv file.
        Should have a 'Temp' column, an 'Rg' column, and and 'Rg_Error' column.
        Rg_Error values can be empty if the error is unknown

    Returns
    ----------
    np.ndarray (2 x m)
        col 0 is temperature and col 1 is Rg for a single chain coil-globule transition.

    """

    if filename[len(filename)-3:] == 'xls':
        our_df = pd.read_excel(filename)

    elif filename[len(filename)-3:] == 'csv':
        our_df = pd.read_csv(filename, sep=r'\s*,\s*', encoding='ascii', engine='python')
    else:
        raise fireballexceptions.FireballException('\n\nERROR: Unable to parse file extension, should be .xls (excel) or .csv')
    
    try:
        data_array = our_df[['Temp','Rg', 'Rg_Error']].to_numpy()
    except Exception as e:
        print('Unable to parse input file')
        print(e)
        raise fireballexceptions.FireballException('Error parsing file')
    
    data_array = data_array[~np.isnan(data_array).transpose()[0]]
    
    return data_array
