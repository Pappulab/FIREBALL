## General configuration

# We define five parameters that will be used as we construct the partial binodal.
# dilute_first_phi and dense_first_phi are used to find the first phi value along
# the dilute and dense arm, respectively.
# dilute_intermediate_increment and dense_intermediate_increment are used to determine
# how many intermediate points need to be found along the dilute and dense arm, respectively.
# partial_binodal_increment defines the allowed distance between independent variable values
# when calculating binodal values.
dilute_first_phi = 0.00001
dense_first_phi = 0.000001
dilute_intermediate_increment = 0.1
dense_intermediate_increment = 0.005
partial_binodal_increment = 10

# The y-variable for our phase diagrams.
y_variable = 'Temperature'
# The label for the y-axis.
y_label = 'Temperature (°C)'

# How to scale the temperature for plotting. If the input is in Kelvin and you want to plot
# in Celsius, this should be 273.15. This is the default value in case no value is passed.
temperature_offset = 273.15

# low_increment (for UCST) and high_increment (for LCST) define the increments for the
# independent variable along the calculated binodal.
low_increment = 1
high_increment = 1

# starting_low_var and starting_high_var define the lowest (for UCST) and highest (for LCST)
# values for the independent variable when calculating the binodal.
starting_low_var = 250
starting_high_var = 350

## Configuration for Gaussian Cluster theory

# Parameters used in Raos and Allegra, 1996
# These include the characteristic ratio, the number of monomers in a Kuhn segment,
# and the length of a single monomer
#char = 10
#ns = 15
#l0 = 3.8

char = 2.5
ns = 7
l0 = .38
