import math
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
from fireball import fireballexceptions, config

'''
Part of the FIREBALL package!

Theory module that defines the Gaussian Cluster theory for
phase separation as described in Raos and Allegra J. Chem Phys 1996.

parameter_list should have the following form: [theta, w2, w3, n]
theta is the theta temperature
w2 is the second virial coefficient of a monomer
w3 is the third virial coefficient of a monomer
n is the polymer length

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

class Theory_Class:
    def __init__(self, params, cur_var):
         
        #Parameter list as defined at the top of this module
        self.params = params
        #The current temperature at which the interpolations are performed
        self.cur_var = cur_var
        #f_scale and g_scale scale the chemical potential (f_scale) and
        #osmotic pressure (g_scale) functions to be on par with one another
        self.f_scale = 1
        self.g_scale = 1000
         
        #Loading in necessary parameters from the config file
        #Characteristic ratio
        self.char = config.char
        #Kuhn length
        self.ns = config.ns
        #Number of kuhn segments
        self.N = self.params[3] / self.ns
        #Second virial coefficient of a kuhn segment
        self.B = self.params[1] * math.sqrt(self.ns)
 
        #The average effective volume per chain monomer in Gaussian Cluster theory.
        #Eq 37 in Raos and Allegra J. Chem Phys 1996
        self.vcbar = math.sqrt(params[2]) * (2 / 3 * math.pi * self.char) ** 1.5
         
        #Eq 41 in Raos and Allegra J. Chem Phys 1996
        self.E = 4 / 3 * math.pi / (2 ** 1.5 * self.ns * self.vcbar)
         
        #Eq 42 in Raos and Allegra J. Chem Phys 1996
        self.F = 4 / 3 * math.pi * self.char ** 1.5 * math.sqrt(self.ns) / (2 ** 1.5 * self.vcbar)
         
        #Create interpolated forms of the necessary free energy functions
        #This is necessary because the compression factor can not be defined explicitly
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
            Temperature (in Kelvin)
     
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
        The real temperature from the reduced temperature in Gaussian Cluster theory.
     
        Parameters
        -------------
        var : float
            Reduced temperature
     
        Returns
        ----------
     
        float
            Returns the temperature in Kelvin
         
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
            Returns sigma
         
        """
        return math.sqrt((1 / math.tanh(q) - 1 / q) / (2 * q))
 
 
# ...........................................................................................
#
#
    def GCT_phi_from_q(self, q):
        """
        The reduced concentration as a function of the compression factor in Gaussian Cluster theory.
        Note that this is different from Eq 20 in Raos and Allegra J. Chem Phys 1996, which is for the single chain
     
        Parameters
        -------------
        q : float
            Compression factor
     
        Returns
        ----------
     
        float
            Returns the reduced concentration
         
        """
        return (-1 * 2 ** 3.5 * self.GCT_sigma2(q) ** 2.5 / 3 * (self.params[2] / (2 * 3 ** 1.5 * self.GCT_sigma2(q) ** 4) + \
                3 * (1 / (2 * math.tanh(q)) + q / (2 * math.sinh(q) ** 2) - 1/q) / (1 / (q * math.sinh(q) ** 2) + \
                1 / (math.tanh(q) * q ** 2) - 2 / (q ** 3))) - self.GCT_T(self.cur_var)) / (self.params[2] * self.F)
 
# ...........................................................................................
#
#
    def GCT_interp_q_from_phi(self):
        """
        The compression factor as a function of the reduced concentration in Gaussian Cluster theory.
      
        Parameters
        -------------
        None
         
        Returns
        ----------
      
        scipy spline class
            Returns the compression factor as a function of phi
          
        """
        
        #y contains the compression factor values
        y = np.logspace(-3, 2, 1000)
        #x contains the reduced phi values from y
        x = [self.GCT_phi_from_q(q) for q in y]
        self.phi_interp_vals = x
        x, y = zip(*sorted(zip(x, y)))
         
        return UnivariateSpline(x, y, k=3, s=0)
  
# ...........................................................................................
#
#
    def GCT_df(self, phi):
        """
        First derivative of the Gaussian Cluster free energy of mixing with respect to
        volume fraction.
     
        Parameters
        -------------
        phi : float
            Volume fraction (must be between 0 and 1)
         
        Returns
        ----------
     
        float
            Returns the instantaneous chemical potential
         
        """
        return (self.GCT_T(self.cur_var) * self.F + self.params[2] * self.F * 2 ** -1.5 * self.GCT_sigma(self.q_interp(self.GCT_phi(phi))) ** -3) * self.GCT_phi(phi) + \
                12 * 3 ** -2.5 * self.params[2] * self.F ** 2 * self.GCT_phi(phi) ** 2 + \
                math.log(self.E * self.GCT_phi(phi)) + self.GCT_T(self.cur_var) * 2 ** -2.5 * self.GCT_sigma(self.q_interp(self.GCT_phi(phi))) ** -3 + self.params[2] / (2 * 3 ** 2.5 * self.GCT_sigma(self.q_interp(self.GCT_phi(phi))) ** 6) + \
                1.5 * (math.log(math.sinh(self.q_interp(self.GCT_phi(phi))) / self.q_interp(self.GCT_phi(phi))) - 0.5 * self.q_interp(self.GCT_phi(phi)) * (1 / math.tanh(self.q_interp(self.GCT_phi(phi))) - 1 / self.q_interp(self.GCT_phi(phi)))) + \
                1.5 * math.log(3 / math.pi) - 3 / 2 - 1.5 * math.log(self.N)
 
# ...........................................................................................
#
#
    def GCT_g(self, phi):
        """
        Osmotic pressure of the Gaussian Cluster free energy of mixing
     
        Parameters
        -------------
        phi : float
            Volume fraction (must be between 0 and 1)
     
        Returns
        ----------
     
        float
            Returns the instantaneous osmotic pressure
         
        """
        return self.GCT_phi(phi) / (self.ns * self.vcbar) + self.F / (2 * self.ns * self.vcbar) * (self.GCT_T(self.cur_var) + \
               self.params[2] * self.GCT_sigma(self.q_interp(self.GCT_phi(phi))) ** -3 / (2 ** 1.5)) * self.GCT_phi(phi) ** 2 + \
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
        self.q_interp = self.GCT_interp_q_from_phi()
        
        #cur_max is necessary to avoid over-extrapolating the splines
        self.cur_max = np.max(self.phi_interp_vals) / np.sqrt(self.N)
        #x contains the reduced phi values used for the extrapolations
        x = np.logspace(-7, np.log10(self.cur_max), 1000)
        df_list = [self.GCT_df(phi) for phi in x]
        g_list = [self.GCT_g(phi) for phi in x]
        
        #df_spline and g_spline are the chemical potential and osmotic pressure, respectively
        self.df_spline = UnivariateSpline(x, df_list, k = 4, s=0)
        self.g_spline = UnivariateSpline(x, g_list, k = 4, s=0)
        
        self.GCT_ddf = self.df_spline.derivative(n=1)
        self.GCT_dg = self.g_spline.derivative(n=1)
      
# ...........................................................................................
#
#
    def build_binodal_function_no_phi(self):
        """
        Function builds the binodal_function, which calculates the square-sum-residuals of the 
        chemical potential (df) and osmotic pressure (g) between two points on the binodal with 
        the same independent variable (temp, salt, etc) values. This version accepts a value
        for the independent variable and leaves the two phi values as unknowns
      
        Parameters
        -------------
        None
         
        Returns
        ----------
      
        1-parameter function
            Returns the binodal_function as a function of x (array-like with form [dense phi value, dilute phi value]) 
      
        """
      
        def f(x):
            return (self.g_scale * self.GCT_g(x[0]) - self.g_scale * self.GCT_g(x[1])) ** 2 + \
                   (self.f_scale * self.GCT_df(x[0]) - self.f_scale * self.GCT_df(x[1])) ** 2
      
        return f

# ...........................................................................................
#
#
    def build_jacobian_function_no_phi(self):
        """
        Function builds the jacobian_function to help determine the binodal points. This incorporates 
        many of the free energy functions above. This version accepts a value for the independent variable
        and leaves the two phi values as unknowns
    
        Parameters
        -------------
        None
        
        Returns
        ----------
    
        1-parameter function
            Returns the jacobian_function as a function of x (array-like with form [dense phi value, dilute phi value]) 
    
        """
    
        def f(x):                                               
            return np.array(((2 * (self.g_scale * self.GCT_g(x[0]) - self.g_scale * self.GCT_g(x[1])) * self.g_scale * self.GCT_dg(x[0])) + \
                            (2 * (self.f_scale * self.GCT_df(x[0]) - self.f_scale * self.GCT_df(x[1])) * self.f_scale * self.GCT_ddf(x[0])), \
                            (2 * (self.g_scale * self.GCT_g(x[0]) - self.g_scale * self.GCT_g(x[1])) * self.g_scale * -1 * self.GCT_dg(x[1])) + \
                            (2 * (self.f_scale * self.GCT_df(x[0]) - self.f_scale * self.GCT_df(x[1]))) * -1 * self.f_scale * self.GCT_ddf(x[1])))
    
        return f

# ...........................................................................................
#
#
    def build_spinodal_function_phi(self):
        """
        Function that builds the spinodal_function as a function of phi. The spinodal is 
        defined by the curve where ddf == 0.
     
        Parameters
        -------------
        None
     
        Returns
        ----------
     
        1-parameter function
            Returns the second derivative of the free energy as a function of phi
     
        """
         
        def f(x):
            return self.f_scale * self.GCT_ddf(x)
         
        return f
       



