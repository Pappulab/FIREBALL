import math
import numpy as np
import pandas as pd
from scipy.interpolate import SmoothBivariateSpline
from fireball import fireballexceptions, config

#from fireball import fireballexceptions, config

'''
Part of the FIREBALL package!

parameter_list should have the following form: [theta, w2, w3, n]
theta is the theta temperature
w2 is the second virial coefficient of a Kuhn segment
w3 is the third virial coefficient of a Kuhn segment
n is the polymer length

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

class Theory_Class:
    def __init__(self, params):
        
        #Parameter list as defined at the top of this module
        self.params = params
        self.f_scale = 1
        self.g_scale = 1000
        
        #Loading in necessary parameters from the config file
        self.char = config.char
        self.ns = config.ns

        self.N = self.params[3] / self.ns
        self.B = self.params[1] * math.sqrt(self.ns)
        self.min_var = config.starting_low_var
        self.max_var = config.starting_high_var

        #The average effective volume per chain monomer in Gaussian Cluster theory.
        #Eq 37 in Raos and Allegra J. Chem Phys 1996
        self.vcbar = math.sqrt(params[2]) * (2 / 3 * math.pi * self.char) ** 1.5
         
        #Eq 41 in Raos and Allegra J. Chem Phys 1996
        self.E = 4 / 3 * math.pi / (2 ** 1.5 * self.ns * self.vcbar)
         
        #Eq 42 in Raos and Allegra J. Chem Phys 1996
        self.F = 4 / 3 * math.pi * self.char ** 1.5 * math.sqrt(self.ns) / (2 ** 1.5 * self.vcbar)
         
        self.Perform_Interpolations()
        
# ...........................................................................................
#
#
    def GCT_phi(self, phi):
        """
        The reduced concentration in Gaussian Cluster theory.
     
        Parameters
        -------------
        phi : float
            Volume fraction (must be between 0 and 1)
     
        Returns
        ----------
     
        float
            Returns the reduced concentration
         
        """
        return phi * self.N ** 0.5
     
# ...........................................................................................
#
#
    def GCT_T(self, var):
        """
        The reduced temperature in Gaussian Cluster theory.
     
        Parameters
        -------------
        var : float
            Temperature
     
        params : list
            list of parameters with form as written at the heading of this module
     
        Returns
        ----------
     
        float
            Returns the reduced temperature
         
        """
        return (var - self.params[0]) / var * self.B * self.N ** 0.5
 

# ...........................................................................................
#
#
    def GCT_T_inverse(self, var):
        """
        The temperature from the reduced temperature in Gaussian Cluster theory.
     
        Parameters
        -------------
        var : float
            Reduced temperature
     
        Returns
        ----------
     
        float
            Returns the real temperature
         
        """
        return self.params[0] * self.B * self.N ** 0.5 / (self.B * self.N ** 0.5 - var)
 
# ...........................................................................................
#
#
    def GCT_sigma2(self, q):
        """
        sigma squared from Gaussian Cluster theory.
        Eq 11 in Raos and Allegra J. Chem Phys 1996
     
        Parameters
        -------------
        q : float
            Compression factor
     
        Returns
        ----------
     
        float
            Returns sigma squared
         
        """
        return (1 / math.tanh(q) - 1 / q) / (2 * q)

# ...........................................................................................
#
#
    def GCT_sigma(self, q):
        """
        sigma from Gaussian Cluster theory.
        Eq 11 in Raos and Allegra J. Chem Phys 1996
    
        Parameters
        -------------
        q : float
            Compression factor
    
        Returns
        ----------
    
        float
            Returns sigma squared
        
        """
        return math.sqrt((1 / math.tanh(q) - 1 / q) / (2 * q))


# ...........................................................................................
#
#
    def GCT_phi_from_q(self, q, var):
        """
        The reduced concentration as a function of the compression factor in Gaussian Cluster theory.
        Note that this is different from Eq 20 in Raos and Allegra J. Chem Phys 1996, which is for the single chain
    
        Parameters
        -------------
        q : float
            Compression factor
    
        var : float
            Temperature
    
        Returns
        ----------
    
        float
            Returns the reduced concentration
        
        """
        return (-1 * 2 ** 3.5 * self.GCT_sigma2(q) ** 2.5 / 3 * (self.params[2] / (2 * 3 ** 1.5 * self.GCT_sigma2(q) ** 4) + \
                3 * (1 / (2 * math.tanh(q)) + q / (2 * math.sinh(q) ** 2) - 1/q) / (1 / (q * math.sinh(q) ** 2) + \
                1 / (math.tanh(q) * q ** 2) - 2 / (q ** 3))) - self.GCT_T(var)) / (self.params[2] * self.F)

# ...........................................................................................
#
#
    def GCT_T_from_q(self, q, phi):
        """
        The reduced temperature as a function of the compression factor in Gaussian Cluster theory.
        Note that this is different from Eq 20 in Raos and Allegra J. Chem Phys 1996, which is for the single chain
    
        Parameters
        -------------
        q : float
            Compression factor
    
        phi : float
            Volume fraction (must be between 0 and 1)
    
        Returns
        ----------
    
        float
            Returns the reduced temperature
        
        """
        return (-1 * 2 ** 3.5 * self.GCT_sigma2(q) ** 2.5 / 3 * (self.params[2] / (2 * 3 ** 1.5 * self.GCT_sigma2(q) ** 4) + \
                3 * (1 / (2 * math.tanh(q)) + q / (2 * math.sinh(q) ** 2) - 1/q) / (1 / (q * math.sinh(q) ** 2) + \
                1 / (math.tanh(q) * q ** 2) - 2 / (q ** 3))) - self.GCT_phi(phi) * self.params[2] * self.F)

# ...........................................................................................
#
#
    def GCT_bivariate_interp_q_from_phi_T(self):
        """
        The compression factor as a function of the reduced temperature in Gaussian Cluster theory.
    
        Parameters
        -------------
        None
    
        Returns
        ----------
        None
        
        """
        
        y = np.linspace(self.min_var, self.max_var, 1000)
        z = np.linspace(.0001, 10, 1000)
        x = [self.GCT_phi_from_q(q, var) for (q, var) in zip(z, y)]
        self.phi_interp_vals = x
        
        return SmoothBivariateSpline(x, y, z, kx=3, ky=3, s=0)
        
# ...........................................................................................
#
#
    def GCT_df(self, phi, var):
        """
        First derivative of the Gaussian Cluster free energy of mixing with respect to
        volume fraction.
    
        Parameters
        -------------
        phi : float
            Volume fraction (must be between 0 and 1)
    
        var : float
            Independent variable (temp, salt, etc.)
        
        Returns
        ----------
    
        float
            Returns the instantaneous chemical potential
        
        """
        return (self.GCT_T(var) * self.F + self.params[2] * self.F * 2 ** -1.5 * self.GCT_sigma(self.q_interp.ev(self.GCT_phi(phi), var)) ** -3) * self.GCT_phi(phi) + \
                12 * 3 ** -2.5 * self.params[2] * self.F ** 2 * self.GCT_phi(phi) ** 2 + \
                math.log(self.E * self.GCT_phi(phi)) + self.GCT_T(var) * 2 ** -2.5 * self.GCT_sigma(self.q_interp.ev(self.GCT_phi(phi), var)) ** -3 + self.params[2] / (2 * 3 ** 2.5 * self.GCT_sigma(self.q_interp.ev(self.GCT_phi(phi), var)) ** 6) + \
                1.5 * (math.log(math.sinh(self.q_interp.ev(self.GCT_phi(phi), var)) / self.q_interp.ev(self.GCT_phi(phi), var)) - 0.5 * self.q_interp.ev(self.GCT_phi(phi), var) * (1 / math.tanh(self.q_interp.ev(self.GCT_phi(phi), var)) - 1 / self.q_interp.ev(self.GCT_phi(phi), var))) + \
                1.5 * math.log(3 / math.pi) - 3 / 2 - 1.5 * math.log(self.N)

# ...........................................................................................
#
#
    def GCT_g(self, phi, var):
        """
        Osmotic pressure of the Gaussian Cluster free energy of mixing
    
        Parameters
        -------------
        phi : float
            Volume fraction (must be between 0 and 1)
    
        var : float
            Independent variable (temp, salt, etc.)
    
        Returns
        ----------
    
        float
            Returns the instantaneous chemical potential
        
        """
        return self.GCT_phi(phi) / (self.ns * self.vcbar) + self.F / (2 * self.ns * self.vcbar) * (self.GCT_T(var) + \
               self.params[2] * self.GCT_sigma(self.q_interp.ev(self.GCT_phi(phi), var)) ** -3 / (2 ** 1.5)) * self.GCT_phi(phi) ** 2 + \
               8 * self.params[2] * self.F ** 2 / (3 ** 2.5 * self.ns * self.vcbar) * self.GCT_phi(phi) ** 3

# ...........................................................................................
#
#
    def Perform_Interpolations(self):
        """
        Performs interpolations to define required functions
        
        Parameters
        -------------
        None
    
        Returns
        ----------
        None
        
        """
        self.q_interp = self.GCT_bivariate_interp_q_from_phi_T()
        
        self.cur_max = np.max(self.phi_interp_vals) / np.sqrt(self.N)
        x = np.logspace(-7, np.log10(self.cur_max), 1000)        
        y = np.linspace(self.min_var, self.max_var, 1000)
        
        df_list = [self.GCT_df(phi, var) for phi, var in zip(x, y)]
        g_list = [self.GCT_g(phi, var) for phi, var in zip(x, y)]
        
        self.df_spline = SmoothBivariateSpline(x, y, df_list, kx=3, ky=3, s=0)
        self.g_spline = SmoothBivariateSpline(x, y, g_list, kx=3, ky=3, s=0)
        
        self.GCT_ddf = self.df_spline.partial_derivative(1, 0)
        self.GCT_dddf = self.df_spline.partial_derivative(2, 0)
        self.GCT_dg = self.g_spline.partial_derivative(1, 0)
    
# ...........................................................................................
#
#
    def build_binodal_function_no_phi(self, var):
        """
        Function builds the binodal_function, which calculates the square-sum-residuals of the 
        chemical potential (df) and osmotic pressure (g) between two points on the binodal with 
        the same independent variable (temp, salt, etc) values. This version accepts a value
        for the independent variable and leaves the two phi values as unknowns
    
        Parameters
        -------------
        var : float
            Independent variable (temp, salt, etc.)
    
        Returns
        ----------
    
        float
            Returns the binodal_function as a function of x (array-like with form [dense phi value, dilute phi value]) 
    
        """
    
        def f(x):
            return (self.g_scale * self.GCT_g.ev(x[0], var) - self.g_scale * self.GCT_g.ev(x[1], var)) ** 2 + \
                   (self.f_scale * self.GCT_df.ev(x[0], var) - self.f_scale * self.GCT_df.ev(x[1], var)) ** 2
    
        return f

# ...........................................................................................
#
#
    def build_jacobian_function_no_phi(self, var):
        """
        Function builds the jacobian_function to help determine the binodal points. This incorporates 
        many of the free energy functions above. This version accepts a value for the independent variable
        and leaves the two phi values as unknowns
     
        Parameters
        -------------
        var : float
            Independent variable (temp, salt, etc.)
     
        Returns
        ----------
     
        float
            Returns the jacobian_function as a function of x (array-like with form [dense phi value, dilute phi value]) 
     
        """
     
        def f(x):                                               
            return np.array(((2 * (self.g_scale * self.GCT_g.ev(x[0], var) - self.g_scale * self.GCT_g.ev(x[1], var)) * self.g_scale * self.GCT_dg(x[0], var)) + \
                            (2 * (self.f_scale * self.GCT_df.ev(x[0], var) - self.f_scale * self.GCT_df.ev(x[1], var)) * self.f_scale * self.GCT_ddf(x[0], var)), \
                            (2 * (self.g_scale * self.GCT_g.ev(x[0], var) - self.g_scale * self.GCT_g.ev(x[1], var)) * self.g_scale * -1 * self.GCT_dg(x[1], var)) + \
                            (2 * (self.f_scale * self.GCT_df.ev(x[0], var) - self.f_scale * self.GCT_df.ev(x[1], var))) * self.f_scale * -1 * self.GCT_ddf(x[1], var)))
     
        return f

# ...........................................................................................
#
#    
    def build_critical_point_function(self):
        """
        Function that builds the critical_point_function. This internal function finds the value of 
        crit_phi and crit_T by solving for when dddf and ddf are zero.    
    
        Parameters
        -------------
        None
        
        Returns
        ----------

        1-parameter function
            Returns the sum of the square second and third derivatives of the potential energy as a function of x
            (array-like with form [phi, temperature])
    
        """

        def f(x):
            return self.f_scale * self.GCT_dddf(x[0], x[1]) ** 2 + self.f_scale * self.GCT_ddf(x[0], x[1]) ** 2
    
        return f

# ...........................................................................................
#
#
    def build_spinodal_function_phi(self, var):
        """
        Function that builds the spinodal_function as a function of phi. The spinodal is 
        defined by the curve where ddf == 0.
    
        Parameters
        -------------
        var : float
            Independent variable (temp, salt, etc.)
    
        Returns
        ----------
    
        1-parameter function
            Returns the second derivative of the potential energy as a function of phi
    
        """
        
        def f(x):
            return self.f_scale * self.GCT_ddf(x, var)
