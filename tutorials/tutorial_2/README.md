FIREBALL Tutorial 2
==============================

## About
This tutorial walks through the process of fitting Gaussian cluster theory to a measured phase diagram of A1-LCD as shown in Figure 2 of the original FIREBALL paper. Data are taken from Bremer et al. (2022).

## Tutorial
From the command-line, `cd` into this directory. Next, enter the following command:

`fireball-fit --filename WT_A1LCD.xls --mode GCT --free_param_list 1_1_1_1 --initial_guess 460_0.0168_0.00632_2000`

Here, `fireball-fit` is the script used for fitting binodals. `--filename WT_A1LCD.xls` tells the program which file the data are in. `--mode GCT` tells the program to use Gaussian cluster theory. `--free_param_list 1_1_1_1` tells the program that `theta`, `w2`, `w3`, and `n` are all free parameters. `--initial_guess ` tells the program how to initialize the values for `theta`, `w2`, `w3`, and `n`.

After running this command, FIREBALL should output a series of "guesses" for the parameters and binodal error values associated with these guesses. Eventually, the optimization should terminate and a plot of the fitted binodal should output to the screen. In addition, the binodal, spinodal, and parameter values should now be saved in this directory.
