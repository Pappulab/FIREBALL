import math
import numpy as np
import pandas as pd
from fireball import fireballexceptions

'''
Part of the FIREBALL package!

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

# ...........................................................................................
#parameter_list should have the following form: [w1, w3, n]
#
def flory_huggins_3B_chi(var, w1):
    """
    Flory-Huggins chi parameter:

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    Returns
    ----------

    float
        Returns the instantaneous Flory-Huggins chi parameter
    
    """
    return w1 / var


def flory_huggins_3B_dchi(var, w1):
    """
    First derivative of the Flory-huggins chi parameter with respect to
    the independent variable (temp, salt, etc.):

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    Returns
    ----------

    float
        Returns the first derivative of the Flory-Huggins chi parameter
    
    """
    return - w1 / (var ** 2)

def flory_huggins_3B_f(phi, var, w1, w3, n):
    """
    Flory-huggins free energy of mixing:

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    Returns
    ----------

    float
        Returns the instantaneous potential energy as defined by FH
    
    """
    return (1 - phi) * math.log(1 - phi) + phi / n * math.log(phi) + flory_huggins_3B_chi(var, w1) * (1 - phi) * phi + w3 * phi ** 3 - (phi ** 3) / 6



# ...........................................................................................
#
#
def flory_huggins_3B_df(phi, var, w1, w3, n):
    """
    First derivative of the Flory-huggins free energy of mixing with respect to
    volume fraction.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    Returns
    ----------

    float
        Returns the first derivative of the potential energy as defined by FH
    
    """

    return math.log(phi) / n + 1 / n - math.log(1 - phi) - 1 + flory_huggins_3B_chi(var, w1) * (1 - 2 * phi) + 3 * w3 * phi ** 2 - (phi ** 2) / 2



# ...........................................................................................
#
#
def flory_huggins_3B_ddf(phi, var, w1, w3, n):
    """
    Second derivative of the Flory-huggins free energy of mixing with respect to
    volume fraction.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    Returns
    ----------

    float
        Returns the second derivative of the potential energy as defined by FH
    
    """

    return 1 / (phi * n) + 1 / (1 - phi) - 2 * flory_huggins_3B_chi(var, w1) + 6 * w3 * phi - phi



# ...........................................................................................
#
#
def flory_huggins_3B_dddf(phi, w3, n):
    """
    Third derivative of the Flory-huggins free energy of mixing with respect to
    volume fraction.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    Returns
    ----------

    float
        Returns the second derivative of the potential energy as defined by FH
    
    """

    return -1 / (n * phi ** 2) + 1 / ((1 - phi) ** 2) + 6 * w3 - 1



# ...........................................................................................
#
#
def flory_huggins_3B_g(phi, var, w1, w3, n):
    """
    Function g finds the y-intercept (free energy) of the free energy tangent line (chemical potential). 
    This is also the osmotic pressure of the system.
    

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    Returns
    ----------

    float
        Returns the y-intercept of the free-energy tangent line 


    """
           
    return math.log(1 - phi) + phi - phi / n + flory_huggins_3B_chi(var, w1) * phi ** 2 - 2 * w3 * phi ** 3 + (phi ** 3) / 3



# ...........................................................................................
#
#
def flory_huggins_3B_dg(phi, var, w1, w3, n):
    """
    Function dg is the derivative of g with respect to phi. This is used to pass the Jacobian 
    to the optimization method while constructing the binodal.

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    Returns
    ----------

    float
        Returns the derivative of the osmotic pressure with respect to phi

    """

    return -1 / (1 - phi) + 1 - 1 / n + 2 * flory_huggins_3B_chi(var, w1) * phi - 6 * w3 * phi ** 2 + phi ** 2 



# ...........................................................................................
#
#
def flory_huggins_3B_df_dvar(phi, var, w1):
    """
    The mixed second derivative of the Flory-huggins free energy of mixing with
    respect to phi and the independent variable (temp, salt, etc.)

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    Returns
    ----------

    float
        Returns the mixed second derivative of the potential energy as defined by FH
    
    """

    return flory_huggins_3B_dchi(var, w1) * (1 - 2 * phi)



# ...........................................................................................
#
#
def flory_huggins_3B_g_dvar(phi, var, w1):
    """
    The mixed second derivative of the g function with
    respect to phi and the independent variable (temp, salt, etc.)

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    Returns
    ----------

    float
        Returns the mixed second derivative of the g function as defined by FH
    
    """

    return flory_huggins_3B_dchi(var, w1) * (phi ** 2)


# ...........................................................................................
#
#
def build_binodal_function(w1, w3, n, phi):
    """
    Function builds the binodal_function, which calculates the square-sum-residuals of the 
    chemical potential (df) and osmotic pressure (g) between two points on the binodal with 
    the same independent variable (temp, salt, etc) values

    Parameters
    -------------
    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    phi : float
        Volume fraction (must be between 0 and 1)

    Returns
    ----------

    float
        Returns the binodal_function as a function of x (array-like with form [phi value, independent variable value]) 

    """

    def f(x):
        return (flory_huggins_3B_g(x[0], x[1], w1, w3, n) - flory_huggins_3B_g(phi, x[1], w1, w3, n)) ** 2 + (flory_huggins_3B_df(x[0], x[1], w1, w3, n) - flory_huggins_3B_df(phi, x[1], w1, w3, n)) ** 2

    return f



# ...........................................................................................
#
#
def build_jacobian_function(w1, w3, n, phi):
    """
    Function builds the jacobian_function to help determine the binodal points. This incorporates 
    many of the free energy functions above.

    Parameters
    -------------
    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    phi : float
        Volume fraction (must be between 0 and 1)

    Returns
    ----------

    float
        Returns the jacobian_function as a function of x (array-like with form [phi value, independent variable value]) 

    """

    def f(x):
        return np.array(((2 * (flory_huggins_3B_g(x[0], x[1], w1, w3, n) - flory_huggins_3B_g(phi, x[1], w1, w3, n)) * flory_huggins_3B_dg(x[0], x[1], w1, w3, n)) + \
                        (2 * (flory_huggins_3B_df(x[0], x[1], w1, w3, n) - flory_huggins_3B_df(phi, x[1], w1, w3, n)) * flory_huggins_3B_ddf(x[0], x[1], w1, w3, n)), \
                        ((2 * (flory_huggins_3B_g(x[0], x[1], w1, w3, n) - flory_huggins_3B_g(phi, x[1], w1, w3, n))) * (flory_huggins_3B_g_dvar(x[0], x[1], w1) - flory_huggins_3B_g_dvar(phi, x[1], w1))) + \
                        ((2 * (flory_huggins_3B_df(x[0], x[1], w1, w3, n) - flory_huggins_3B_df(phi, x[1], w1, w3, n))) * (flory_huggins_3B_df_dvar(x[0], x[1], w1) - flory_huggins_3B_df_dvar(phi, x[1], w1)))))

    return f


# ...........................................................................................
#
#
def build_binodal_function_no_phi(w1, w3, n, var):
    """
    Function builds the binodal_function, which calculates the square-sum-residuals of the 
    chemical potential (df) and osmotic pressure (g) between two points on the binodal with 
    the same independent variable (temp, salt, etc) values. This version accepts a value
    for the independent variable and leaves the two phi values as unknowns

    Parameters
    -------------
    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    var : float
        Independent variable (temp, salt, etc.)

    Returns
    ----------

    float
        Returns the binodal_function as a function of x (array-like with form [dense phi value, dilute phi value]) 

    """

    def f(x):
        return (flory_huggins_3B_g(x[0], var, w1, w3, n) - flory_huggins_3B_g(x[1], var, w1, w3, n)) ** 2 + (flory_huggins_3B_df(x[0], var, w1, w3, n) - flory_huggins_3B_df(x[1], var, w1, w3, n)) ** 2

    return f


# ...........................................................................................
#
#
def build_jacobian_function_no_phi(w1, w3, n, var):
    """
    Function builds the jacobian_function to help determine the binodal points. This incorporates 
    many of the free energy functions above. This version accepts a value for the independent variable
    and leaves the two phi values as unknowns

    Parameters
    -------------
    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    var : float
        Independent variable (temp, salt, etc.)

    Returns
    ----------

    float
        Returns the jacobian_function as a function of x (array-like with form [dense phi value, dilute phi value]) 

    """

    def f(x):                                               
        return np.array(((2 * (flory_huggins_3B_g(x[0], var, w1, w3, n) - flory_huggins_3B_g(x[1], var, w1, w3, n)) * flory_huggins_3B_dg(x[0], var, w1, w3, n)) + \
                        (2 * (flory_huggins_3B_df(x[0], var, w1, w3, n) - flory_huggins_3B_df(x[1], var, w1, w3, n)) * flory_huggins_3B_ddf(x[0], var, w1, w3, n)), \
                        (2 * (flory_huggins_3B_g(x[0], var, w1, w3, n) - flory_huggins_3B_g(x[1], var, w1, w3, n)) * -flory_huggins_3B_dg(x[1], var, w1, w3, n)) + \
                        (2 * (flory_huggins_3B_df(x[0], var, w1, w3, n) - flory_huggins_3B_df(x[1], var, w1, w3, n))) * -flory_huggins_3B_ddf(x[1], var, w1, w3, n)))

    return f


# ...........................................................................................
#
#    
def build_critical_phi_function(w3, n):
    """
    Function that builds the critical_phi_function. This internal function finds the value of 
    crit_phi by solving for when dddf is zero. This value is independent of the independent variable.    

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    w3 : float
        Three-body interaction term

    Returns
    ----------

    1-parameter function
        Returns the third derivative of the potential energy as a function of phi
    

    """

    def f(x):
        return flory_huggins_3B_dddf(x[0], w3, n)

    return f

# ...........................................................................................
#
#
def build_spinodal_function_phi(var, w1, w3, n):
    """
    Function that builds the spinodal_function as a function of phi. The spinodal is 
    defined by the curve where ddf == 0. If crit_phi is known, this equation can be 
    solved to find the critical point.    

    Parameters
    -------------
    var : float
        Independent variable (temp, salt, etc.)

    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    Returns
    ----------

    1-parameter function
        Returns the second derivative of the potential energy as a function of phi
    

    """

    def f(x):
        return flory_huggins_3B_ddf(x, var, w1, w3, n)

    return f

# ...........................................................................................
#
#
def build_spinodal_function_var(phi, w1, w3, n):
    """
    Function that builds the spinodal_function as a function of the independent variable 
    (temp, salt, etc). The spinodal is defined by the curve where ddf == 0. If crit_phi 
    is known, this equation can be solved to find the critical point.    

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization

    Returns
    ----------

    1-parameter function
        Returns the second derivative of the potential energy as a function of the independent variable (temp, salt, etc)
    

    """

    def f(x):
        return flory_huggins_3B_ddf(phi, x, w1, w3, n)

    return f


# ...........................................................................................
#
#
def calculate_spinodal_var(phi, w1, w3, n):
    """
    Function that finds the value of the independent variable (temp, salt, etc) on the spinodal 
    for a given phi value by explicitly solving for when ddf is zero. Substituting crit_phi for phi
    will return the critical value for the independent variable 

    Parameters
    -------------
    phi : float
        Volume fraction (must be between 0 and 1)

    w1 : float
        Flory chi parameter energy term

    w3 : float
        Three-body interaction term

    n : int
        Degree of polymerization
        
    Returns
    ----------

    float
        Returns the value of the independent variable (temp, salt, etc) on the spinodal for the given phi value
    

    """

    return (2 * w1) / (1 / (phi * n) + 1 / (1 - phi) + 6 * w3 * phi - phi)

