FIREBALL Tutorial 1
==============================

## About
This tutorial walks through the process of fitting Flory-Huggins theory to a measured phase diagram of A1-LCD as shown in Figure 1 of the original FIREBALL paper. Data are taken from Bremer et al. (2022).

## Tutorial
From the command-line, `cd` into this directory. Next, enter the following command:

`fireball-fit --filename WT_A1LCD.xls --free_param_list 1_1_0 --initial_guess 220_.4_137`

Here, `fireball-fit` is the script used for fitting binodals. `--filename WT_A1LCD.xls` tells the program which file the data are in. `--free_param_list 1_1_0` tells the program that `u` and `w3` are free parameters, whereas `n` is fixed. `--initial_guess 220_.4_137` tells the program how to initialize the values for `u`, `w3`, and `n`.

After running this command, FIREBALL should output a series of "guesses" for the parameters and binodal error values associated with these guesses. Eventually, the optimization should terminate and a plot of the fitted binodal should output to the screen. In addition, the binodal, spinodal, and parameter values should now be saved in this directory.