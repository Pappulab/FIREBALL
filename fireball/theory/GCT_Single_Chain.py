import math
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
from fireball import fireballexceptions, config

'''
Part of the FIREBALL package!

Theory module that defines the Gaussian Cluster Theory for
single chains as described in Raos and Allegra J. Chem Phys 1996.

parameter_list should have the following form: [theta, w2, w3, n]
theta is the theta temperature
w2 is the second virial coefficient of a Kuhn segment
w3 is the third virial coefficient of a Kuhn segment
n is the degree of polymerization

Originally created on Sep 30, 2019 by Mina Farag
Updated Spring 2021 by Alex Holehouse and Mina Farag

'''

class Theory_Class:
    def __init__(self, params, partial_check = False):
         
        #Parameter list as defined at the top of this module
        self.params = params
         
        #Loading in necessary parameters from the config file
        
        #Characteristic ratio
        self.char = config.char
        #Kuhn length
        self.ns = config.ns
        #Length of a monomer
        self.l0 = config.l0
        #Number of kuhn segments
        self.N = self.params[3] / self.ns
        #Second virial coefficient of a kuhn segment
        self.B = self.params[1] * math.sqrt(self.ns)
        #Describes the size of a kuhn segment
        self.l2 = self.ns * self.char * self.l0 ** 2
 
        #The average effective volume per chain monomer in Gaussian Cluster theory.
        #Eq 37 in Raos and Allegra J. Chem Phys 1996
        self.vcbar = math.sqrt(params[2]) * (2 / 3 * math.pi * self.char) ** 1.5
         
        #Eq 41 in Raos and Allegra J. Chem Phys 1996
        self.E = 4 / 3 * math.pi / (2 ** 1.5 * self.ns * self.vcbar)
         
        #Eq 42 in Raos and Allegra J. Chem Phys 1996
        self.F = 4 / 3 * math.pi * self.char ** 1.5 * math.sqrt(self.ns) / (2 ** 1.5 * self.vcbar)
        
        #Create interpolated forms of the necessary free energy functions
        #This is necessary because the compression factor can not be defined explicitly
        #If partial_check is set to False, the full coil-globule transition is calculated
        self.Perform_Interpolations(partial_check)
        
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
    def GCT_rg2(self, q):
        """
        <Rg^2> from Gaussian Cluster theory.
        Eq 11 in Raos and Allegra J. Chem Phys 1996
     
        Parameters
        -------------
        q : float
            Compression factor
     
        Returns
        ----------
     
        float
            Returns <rg^2>
         
        """
        return self.GCT_sigma2(q) * self.N * self.l2

# ...........................................................................................
#
#
    def GCT_alpha(self, rg2):
        """
        alpha from Gaussian Cluster theory.
        Eq 26 in Raos and Allegra J. Chem Phys 1996
     
        Parameters
        -------------
        rg2 : float
            <Rg^2>
     
        Returns
        ----------
     
        float
            Returns alpha
         
        """
        return np.sqrt(6 * rg2 / (self.N * self.l2))

# ...........................................................................................
#
#
    def GCT_T_from_q(self, q):
        """
        The reduced temperature as a function of the single-chain compression factor in Gaussian Cluster theory.
        This is the same as Eq 20 in Raos and Allegra J. Chem Phys 1996
     
        Parameters
        -------------
        q : float
            Compression factor
     
        Returns
        ----------
     
        float
            Returns the reduced temperature
         
        """
        return -1 * 2 ** 3.5 * self.GCT_sigma2(q) ** 2.5 / 3 * (self.params[2] / (2 * 3 ** 1.5 * self.GCT_sigma2(q) ** 4) + \
               3 * (1 / (2 * math.tanh(q)) + q / (2 * math.sinh(q) ** 2) - 1/q) / (1 / (q * math.sinh(q) ** 2) + \
               1 / (math.tanh(q) * q ** 2) - 2 / (q ** 3)))

# ...........................................................................................
#
#
    def GCT_interp_q_from_T(self):
        """
        The compression factor as a function of the reduced temperature in Gaussian Cluster theory.
      
        Parameters
        -------------
        None
         
        Returns
        ----------
      
        scipy spline class
            Returns the compression factor as a function of the reduced temperature
          
        """
        
        #y contains the compression factor values
        y = np.logspace(-3.5, 2, 1000)
        #x contains the reduced temperature values from y
        x = [self.GCT_T_from_q(q) for q in y]
        self.T_interp_vals = x
        x, y = zip(*sorted(zip(x, y)))
         
        return UnivariateSpline(x, y, k=3, s=0)
  
 
# ...........................................................................................
#
#
    def Perform_Interpolations(self, partial_check):
        """
        Performs interpolations to define required functions
         
        Parameters
        -------------
        partial_check : boolean
            Flag which, if set to True, stops the program from calculating
            the full coil-globule transition
     
        Returns
        ----------
        None
         
        """    
        self.q_interp = self.GCT_interp_q_from_T()
        
        #cur_min and cur_max are necessary to avoid over-extrapolating the spline 
        self.cur_min = np.min(self.T_interp_vals)
        self.cur_max = np.max(self.T_interp_vals)
        
        #If partial_check is set to False, the full coil-globule transition is calculated
        if not partial_check:
            self.reduced_t_vals = np.linspace(self.cur_min, self.cur_max, 1000)
            self.real_t_vals = [self.GCT_T_inverse(var) for var in self.reduced_t_vals]
            self.rg_vals = [self.GCT_rg2(self.q_interp(var)) for var in self.reduced_t_vals]
            self.alpha_vals = [self.GCT_alpha(rg2) for rg2 in self.rg_vals]
      
       



