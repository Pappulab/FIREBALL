#!/usr/bin/env python

###############################################

import argparse

import fireball
from fireball import full_single_chain, io_single_chain

def startup():
    print('Running FIREBALL-SINGLE-CHAIN-DRAW')


## ===================================================================================================
##                              Main Script - hold onto your hat!
## ===================================================================================================

def main():
    # this builds the parser object, which we'll use to add in command-line options
    parser = argparse.ArgumentParser()

    # this is where command-line options are defined
    parser.add_argument("--filename", help="[OPTIONAL] Input file with information for plotting data")
    parser.add_argument("--param_list", help="[OPTIONAL] parameter values (format X_Y_Z where X,Y,Z depend on the fitting procedure. Run --fitting-info for more) [Default = 345_.02_.0076_3000]",  default = "345_.02_.0076_3000")
    parser.add_argument("--temperature-offset", help="[OPTIONAL] Value that defines the temperature offset for plotting [default = 273.15, assumes data should be plotted in Celsius]", default='273.15')
    parser.add_argument("--LCST-check", help="[OPTIONAL] Flag which, if set to true, tells the algorithm that this is a lower critical system (default = False)", action='store_true')
    parser.add_argument("--increment", help="[OPTIONAL] The temperature increment [default=1.0]", default='1.0')
    parser.add_argument("--mode", help="[OPTIONAL] Defines the fitting mode to be used. Valid options are: GCT_Single_Chain [default=GCT_Single_Chain]", default='GCT_Single_Chain')
    parser.add_argument("--outname", help="[OPTIONAL] Name of output file where single chain data are written", default='single-chain.csv')
    parser.add_argument("--noplot", help="Flag which, if provided, means no figure is generated to screen but the single chain data are just saved to disk", action='store_true')
    parser.add_argument("--version", help="Flag which, if set to true, prints software version and exits the program", action='store_true')

    args = parser.parse_args()

    if args.version:
        print("Version: %s"%(str(fireball.__version__)))
        exit(0)

    # Mode specific stuff in terms of input options can happen here
    if args.mode == 'GCT_Single_Chain':

        parameter_list = [float(x) for x in args.param_list.split('_')]

        if len(parameter_list) != 4:
            print("parameter_list must have 4 values for the given theory!")

    else:
        print("No valid --mode offered. Valid options are: 'GCT_Single_Chain'")
        exit(1)

    # set some defaults and read in data if provided
    try:
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
        data_array = io_single_chain.read_file(args.filename)
    else:
        data_array = None

    if args.noplot:
        plotter = False
    else:
        plotter = True

    full_single_chain.coil_globule_maker(parameter_list, mode, data_array, temperature_offset=temperature_offset, increment=increment, outname=outname, plotter=plotter)

