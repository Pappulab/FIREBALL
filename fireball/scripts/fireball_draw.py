#!/usr/bin/env python

###############################################

import os
import warnings
import argparse
import numpy as np

import fireball
from fireball import full_binodal, io

def startup():
    print('Running FIREBALL-DRAW')


## ===================================================================================================
##                              Main Script - hold onto your hat!
## ===================================================================================================

def main():
    # this builds the parser object, which we'll use to add in command-line options
    parser = argparse.ArgumentParser()

    # this is where command-line options are defined
    parser.add_argument("--filename", help="[OPTIONAL] Input file with information for plotting data")
    parser.add_argument("--param_list", help="[OPTIONAL] parameter values (format X_Y_Z where X,Y,Z depend on the fitting procedure. Run --fitting-info for more) [Default = 209_0.31_100]",  default = "209_0.31_100")
    parser.add_argument("--conversion-factor", help="[OPTIONAL] Conversion factor that converts concentration in mg/ml to volume fraction (default = 1310 mg/ml)", default="1310")
    parser.add_argument("--temperature-offset", help="[OPTIONAL] Value that defines the temperature offset for plotting [default = 273.15, assumes data should be plotted in Celsius]", default='273.15')
    parser.add_argument("--LCST-check", help="[OPTIONAL] Flag which, if set to true, tells the algorithm that this is a lower critical system (default = False)", action='store_true')
    parser.add_argument("--increment", help="[OPTIONAL] The phase diagram increment [default=1.0]", default='1.0')
    parser.add_argument("--mode", help="[OPTIONAL] Defines the fitting mode to be used. Valid options are: 'flory_huggins_3B' and 'GCT' [default=flory_huggins_3B]", default='flory_huggins_3B')
    parser.add_argument("--outname", help="[OPTIONAL] Prefix of output files where binodal/spinodal data are written [default = Draw]", default='Draw')
    parser.add_argument("--noplot", help="Flag which, if provided, means no figure is generated to screen but the binodal data are just saved to disk", action='store_true')
    parser.add_argument("--version", help="Flag which, if set to true, prints software version and exits the program", action='store_true')

    args = parser.parse_args()

    if args.version:
        print("Version: %s"%(str(fireball.__version__)))
        exit(0)

    # Mode specific stuff in terms of input options can happen here
    if args.mode == 'flory_huggins_3B' or args.mode == 'flory_huggins_3B_LCST':

        parameter_list = [float(x) for x in args.param_list.split('_')]

        if len(parameter_list) != 3:
            print("parameter_list must have 3 values for the given theory!")

    elif args.mode == 'GCT' or args.mode == 'GCT_Bivariate':

        parameter_list = [float(x) for x in args.param_list.split('_')]

        if len(parameter_list) != 4:
            print("parameter_list must have 4 values for the given theory!")

    else:
        print("No valid --mode offered. Valid options are: 'flory_huggins_3B' or 'GCT'")
        exit(1)

    # set some defaults and read in data if provided
    try:
        conv = float(args.conversion_factor)
        increment = float(args.increment)
        outname = args.outname
        temperature_offset = float(args.temperature_offset)
        mode = args.mode
    except Exception as e:
        print('\n\n\nSomething went wrong when parsing commandline arguments. See below...')
        print(e)
        exit(1)

    if args.LCST_check:
        regime = 'high'
    else:
        regime = 'low'

    if args.filename:
        (dilute_array, dense_array) = io.read_file(args.filename, conv)
    else:
        dilute_array = None
        dense_array = None

    if args.noplot:
        plotter = False
    else:
        plotter = True

    full_binodal.binodal_maker(parameter_list, mode, regime, dilute_array, dense_array, temperature_offset=temperature_offset, increment=increment, outname=outname, plotter=plotter)


if __name__=="__main__":
    main()
