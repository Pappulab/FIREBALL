import math
import numpy as np
import pandas as pd
from fireball import fireballexceptions

'''
Part of the FIREBALL package!

parameter_list should have the following form: [w1, w3, n]
w1 is the Flory chi parameter energy term
w3 is a three-body interaction term
n is the degree of polymerization

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

# ...........................................................................................
#
#
def flory_huggins_3B_chi(var, params):
    """
    Flory-Huggins chi parameter:

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the instantaneous Flory-Huggins chi parameter
    
    """
    return params[0] * var


def flory_huggins_3B_dchi(var, params):
    """
    First derivative of the Flory-huggins chi parameter with respect to
    the independent variable (temp, salt, etc.):

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the first derivative of the Flory-Huggins chi parameter
    
    """
    return params[0]

def flory_huggins_3B_f(phi, var, params):
    """
    Flory-huggins free energy of mixing:

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the instantaneous potential energy as defined by FH
    
    """
    return (1 - phi) * math.log(1 - phi) + phi / params[2] * math.log(phi) + flory_huggins_3B_chi(var, params) * (1 - phi) * phi + params[1] * phi ** 3 - (phi ** 3) / 6



# ...........................................................................................
#
#
def flory_huggins_3B_df(phi, var, params):
    """
    First derivative of the Flory-huggins free energy of mixing with respect to
    volume fraction.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the first derivative of the potential energy as defined by FH
    
    """

    return math.log(phi) / params[2] + 1 / params[2] - math.log(1 - phi) - 1 + flory_huggins_3B_chi(var, params) * (1 - 2 * phi) + 3 * params[1] * phi ** 2 - (phi ** 2) / 2



# ...........................................................................................
#
#
def flory_huggins_3B_ddf(phi, var, params):
    """
    Second derivative of the Flory-huggins free energy of mixing with respect to
    volume fraction.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the second derivative of the potential energy as defined by FH
    
    """

    return 1 / (phi * params[2]) + 1 / (1 - phi) - 2 * flory_huggins_3B_chi(var, params) + 6 * params[1] * phi - phi



# ...........................................................................................
#
#
def flory_huggins_3B_dddf(phi, params):
    """
    Third derivative of the Flory-huggins free energy of mixing with respect to
    volume fraction.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the second derivative of the potential energy as defined by FH
    
    """

    return -1 / (params[2] * phi ** 2) + 1 / ((1 - phi) ** 2) + 6 * params[1] - 1



# ...........................................................................................
#
#
def flory_huggins_3B_g(phi, var, params):
    """
    Function g finds the y-intercept (free energy) of the free energy tangent line (chemical potential). 
    This is also the osmotic pressure of the system.
    

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the y-intercept of the free-energy tangent line 


    """
           
    return math.log(1 - phi) + phi - phi / params[2] + flory_huggins_3B_chi(var, params) * phi ** 2 - 2 * params[1] * phi ** 3 + (phi ** 3) / 3



# ...........................................................................................
#
#
def flory_huggins_3B_dg(phi, var, params):
    """
    Function dg is the derivative of g with respect to phi. This is used to pass the Jacobian 
    to the optimization method while constructing the binodal.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the derivative of the osmotic pressure with respect to phi

    """

    return -1 / (1 - phi) + 1 - 1 / params[2] + 2 * flory_huggins_3B_chi(var, params) * phi - 6 * params[1] * phi ** 2 + phi ** 2 



# ...........................................................................................
#
#
def flory_huggins_3B_df_dvar(phi, var, params):
    """
    The mixed second derivative of the Flory-huggins free energy of mixing with
    respect to phi and the independent variable (temp, salt, etc.)

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the mixed second derivative of the potential energy as defined by FH
    
    """

    return flory_huggins_3B_dchi(var, params) * (1 - 2 * phi)



# ...........................................................................................
#
#
def flory_huggins_3B_g_dvar(phi, var, params):
    """
    The mixed second derivative of the g function with
    respect to phi and the independent variable (temp, salt, etc.)

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the mixed second derivative of the g function as defined by FH
    
    """

    return flory_huggins_3B_dchi(var, params) * (phi ** 2)


# ...........................................................................................
#
#
def build_binodal_function(phi, params):
    """
    Function builds the binodal_function, which calculates the square-sum-residuals of the 
    chemical potential (df) and osmotic pressure (g) between two points on the binodal with 
    the same independent variable (temp, salt, etc) values

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the binodal_function as a function of x (array-like with form [phi value, independent variable value]) 

    """

    def f(x):
        return (flory_huggins_3B_g(x[0], x[1], params) - flory_huggins_3B_g(phi, x[1], params)) ** 2 + (flory_huggins_3B_df(x[0], x[1], params) - flory_huggins_3B_df(phi, x[1], params)) ** 2

    return f



# ...........................................................................................
#
#
def build_jacobian_function(phi, params):
    """
    Function builds the jacobian_function to help determine the binodal points. This incorporates 
    many of the free energy functions above.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the jacobian_function as a function of x (array-like with form [phi value, independent variable value]) 

    """

    def f(x):
        return np.array(((2 * (flory_huggins_3B_g(x[0], x[1], params) - flory_huggins_3B_g(phi, x[1], params)) * flory_huggins_3B_dg(x[0], x[1], params)) + \
                        (2 * (flory_huggins_3B_df(x[0], x[1], params) - flory_huggins_3B_df(phi, x[1], params)) * flory_huggins_3B_ddf(x[0], x[1], params)), \
                        ((2 * (flory_huggins_3B_g(x[0], x[1], params) - flory_huggins_3B_g(phi, x[1], params))) * (flory_huggins_3B_g_dvar(x[0], x[1], params) - flory_huggins_3B_g_dvar(phi, x[1], params))) + \
                        ((2 * (flory_huggins_3B_df(x[0], x[1], params) - flory_huggins_3B_df(phi, x[1], params))) * (flory_huggins_3B_df_dvar(x[0], x[1], params) - flory_huggins_3B_df_dvar(phi, x[1], params)))))

    return f


# ...........................................................................................
#
#
def build_binodal_function_no_phi(var, params):
    """
    Function builds the binodal_function, which calculates the square-sum-residuals of the 
    chemical potential (df) and osmotic pressure (g) between two points on the binodal with 
    the same independent variable (temp, salt, etc) values. This version accepts a value
    for the independent variable and leaves the two phi values as unknowns

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the binodal_function as a function of x (array-like with form [dense phi value, dilute phi value]) 

    """

    def f(x):
        return (flory_huggins_3B_g(x[0], var, params) - flory_huggins_3B_g(x[1], var, params)) ** 2 + (flory_huggins_3B_df(x[0], var, params) - flory_huggins_3B_df(x[1], var, params)) ** 2

    return f


# ...........................................................................................
#
#
def build_jacobian_function_no_phi(var, params):
    """
    Function builds the jacobian_function to help determine the binodal points. This incorporates 
    many of the free energy functions above. This version accepts a value for the independent variable
    and leaves the two phi values as unknowns

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    float
        Returns the jacobian_function as a function of x (array-like with form [dense phi value, dilute phi value]) 

    """

    def f(x):                                               
        return np.array(((2 * (flory_huggins_3B_g(x[0], var, params) - flory_huggins_3B_g(x[1], var, params)) * flory_huggins_3B_dg(x[0], var, params)) + \
                        (2 * (flory_huggins_3B_df(x[0], var, params) - flory_huggins_3B_df(x[1], var, params)) * flory_huggins_3B_ddf(x[0], var, params)), \
                        (2 * (flory_huggins_3B_g(x[0], var, params) - flory_huggins_3B_g(x[1], var, params)) * -flory_huggins_3B_dg(x[1], var, params)) + \
                        (2 * (flory_huggins_3B_df(x[0], var, params) - flory_huggins_3B_df(x[1], var, params))) * -flory_huggins_3B_ddf(x[1], var, params)))

    return f


# ...........................................................................................
#
#    
def build_critical_phi_function(params):
    """
    Function that builds the critical_phi_function. This internal function finds the value of 
    crit_phi by solving for when dddf is zero. This value is independent of the independent variable.    

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    1-parameter function
        Returns the third derivative of the potential energy as a function of phi
    

    """

    def f(x):
        return flory_huggins_3B_dddf(x[0], params)

    return f

# ...........................................................................................
#
#
def build_spinodal_function_phi(var, params):
    """
    Function that builds the spinodal_function as a function of phi. The spinodal is 
    defined by the curve where ddf == 0. If crit_phi is known, this equation can be 
    solved to find the critical point.    

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    1-parameter function
        Returns the second derivative of the potential energy as a function of phi
    

    """

    def f(x):
        return flory_huggins_3B_ddf(x, var, params)

    return f

# ...........................................................................................
#
#
def build_spinodal_function_var(phi, params):
    """
    Function that builds the spinodal_function as a function of the independent variable 
    (temp, salt, etc). The spinodal is defined by the curve where ddf == 0. If crit_phi 
    is known, this equation can be solved to find the critical point.    

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    params : list
        list of parameters with form as written at the heading of this module

    Returns
    ----------

    1-parameter function
        Returns the second derivative of the potential energy as a function of the independent variable (temp, salt, etc)
    

    """

    def f(x):
        return flory_huggins_3B_ddf(phi, x, params)

    return f

