FIREBALL Tutorial 3
==============================

## About
This tutorial walks through the process of fitting Gaussian cluster theory to a computationally derived coil-to-globule transition of the synthetic protein (QGQSPYG)9 as shown in Figure 3 of the original FIREBALL paper. Data are taken from Zeng et al. (2020).

## Tutorial
From the command-line, `cd` into this directory. Next, enter the following command:

`fireball-single-chain-fit --filename QGQSPYG9_Rg.csv --initial_guess 345_.58_.018_90`

Here, `fireball-single-chain-fit` is the script used for fitting coil-to-globule transition data. `--filename QGQSPYG9_Rg.csv` tells the program which file the data are in. `--initial_guess 345_.58_.018_90` tells the program how to initialize the values for `theta`, `w2`, `w3`, and `n`, respectively.

After running this command, FIREBALL should output a series of "guesses" for the parameters and fitting error values associated with these guesses. Eventually, the optimization should terminate and a plot of the fitted transition should output to the screen. In addition, the fitted transition and parameter values should now be saved in this directory.