#!/usr/bin/env python

###############################################

import os
import warnings
import argparse
import numpy as np

import fireball
from fireball import optimizer_single_chain, full_single_chain, io_single_chain

def startup():
    print('Running FIREBALL-SINGLE-CHAIN-FIT')


## ===================================================================================================
##                              Main Script - hold onto your hat!
## ===================================================================================================

def main():
    # this builds the parser object, which we'll use to add in command-line options
    parser = argparse.ArgumentParser()

    # this is where command-line options are defined
    parser.add_argument("--filename", help="Input file with data")
    parser.add_argument("--initial_guess", help="[OPTIONAL] Initial guess for values (format X_Y_Z where X,Y,Z depend on the fitting procedure) [Default = 345_.02_.0076_3000]",  default = "345_.02_.0076_3000")
    parser.add_argument("--free_param_list", help="[OPTIONAL] Define which parameters should be free vs fixed (format X_Y_Z where X,Y,Z are 0 if the parameter is fixed and 1 if it is free) [default = 1_1_1_1]", default = "1_1_1_1")
    parser.add_argument("--temperature-offset", help="[OPTIONAL] Value that defines the temperature offset for plotting [default = 273.15, assumes data should be plotted in Celsius]", default='273.15')
    parser.add_argument("--LCST-check", help="[OPTIONAL] Flag which, if set to true, tells the algorithm that this is a lower critical system (default = False)", action='store_true')
    parser.add_argument("--fitting-method", help="[OPTIONAL] Method used to fit [default = 'Nelder-Mead']", default="Nelder-Mead")
    parser.add_argument("--maxfev", help="[OPTIONAL] Max number of function evaluations allowed during optimization [default=200]", default='2000')
    parser.add_argument("--xatol", help="[OPTIONAL] The absolute error in the parameter set between iterations that is acceptable for convergence [default=0.001]", default='0.001')
    parser.add_argument("--fatol", help="[OPTIONAL] The absolute error in the square-residuals between iterations that is acceptable for convergence [default=0.001]", default='0.001')
    parser.add_argument("--increment", help="[OPTIONAL] The temperature increment [default=1.0]", default='1.0')
    parser.add_argument("--mode", help="[OPTIONAL] Defines the fitting mode to be used. Valid options are: GCT_Single_Chain [default=GCT_Single_Chain]", default='GCT_Single_Chain')
    parser.add_argument("--outname", help="[OPTIONAL] Name of output file where single chain data are written", default='single-chain.csv')
    parser.add_argument("--outname-fitfile", help="[OPTIONAL] Name of output file where final fitting parameters are written [default = fit_params.csv]", default='fit_params.csv')
    parser.add_argument("--show-partial-fit", help="Flag which, if set to true, will show the coil to globule transition generated from the partial fit", action='store_true')
    parser.add_argument("--noplot", help="Flag which, if provided, means no figure is generated to screen but the single chain data are just saved to disk", action='store_true')
    parser.add_argument("--version", help="Flag which, if set to true, prints software version and exits the program", action='store_true')

    args = parser.parse_args()

    if args.version:
        print("Version: %s"%(str(fireball.__version__)))
        exit(0)

    # Mode specific stuff in terms of input options can happen here
    if args.mode == 'GCT_Single_Chain':

        init_guess = [float(x) for x in args.initial_guess.split('_')]
        free_param_checklist = [int(x) for x in args.free_param_list.split('_')]

        if len(init_guess) != 4:
            print("init_guess must have 4 values for the given theory!")
        if len(free_param_checklist) != 4:
            print("free_param_checklist must have 4 values for the given theory!")

    else:
        print("No valid --mode offered. Valid options are: 'GCT_Single_Chain'")
        exit(1)

    if args.filename is None:
        print('Must pass a filename with --filename flag')
        exit(1)

    ## extract out parameters as passed in
    # ..........................................

    try:
        fitting_method = args.fitting_method # The fitting method to be used for the optimization
        maxfev = int(args.maxfev) # The maximum number of function evaluations allowed during the optimization
        xatol = float(args.xatol) # The absolute error in the parameter set between iterations that is acceptable for convergence
        fatol = float(args.fatol) # The absolute error in the sum of the square-residuals between iterations that is acceptable for convergence
        increment = float(args.increment)
        mode = args.mode
        outname_fitfile = args.outname_fitfile
        outname = args.outname
        temperature_offset = float(args.temperature_offset)
    except Exception as e:
        print('\n\n\nSomething went wrong when parsing commandline arguments. See below...')
        print(e)
        exit(1)

    if args.LCST_check:
        regime = 'high'
    else:
        regime = 'low'

    if args.show_partial_fit:
        partial_fit = True
    else:
        partial_fit = False

    if args.noplot:
        plotter = False
    else:
        plotter = True

    # Read in single chain data
    data_array = io_single_chain.read_file(args.filename)

    # Model-specific fitting routines can be called here
    if args.mode == 'GCT_Single_Chain':

        # stage 1 is fitting to the data
        fit_obj = optimizer_single_chain.optimizer(mode, data_array, init_guess, free_param_checklist, temperature_offset, fitting_method, maxfev, xatol, fatol, partial_fit)

        # stage 2 is plotting the data
        full_parameter_list = []
        counter_free = 0
        for index, element in enumerate(free_param_checklist):
            if element == 1:
                full_parameter_list.append(fit_obj.x[counter_free])
                counter_free += 1
            else:
                full_parameter_list.append(init_guess[index])

        if args.mode == 'GCT_Single_Chain':
            theta = full_parameter_list[0]
            v2 = full_parameter_list[1]
            v3 = full_parameter_list[2]
            n = full_parameter_list[3]

            # write out fitting info
            with open(outname_fitfile,'a+') as fp:
                fp.write('name, theta, v2, v3, n, residuals\n')
                fp.write('%s, %4.4f, %4.4f, %4.4f, %4.4f, %4.4f\n'%(outname, theta, v2, v3, n, fit_obj.fun))

        # finally generate full coil to globule transition!
        full_single_chain.coil_globule_maker(full_parameter_list, mode, data_array, temperature_offset=temperature_offset, increment=increment, outname=outname, plotter=plotter)

    else:
        print('Unknown --mode flag passed')
        exit(1)


if __name__=="__main__":
    main()
