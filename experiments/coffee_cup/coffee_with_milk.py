"""Buoyancy experiment: mix hot coffee with cold milk."""

# Third party imports
import numpy as np

# Local Nyles imports
from nyles import Nyles
from parameters import UserParameters


# Densities in kg m^-3 (the unit does not matter, only the density ratio does)
rho_coffee = 970  # at 85 °C (density of fresh water)
rho_milk = 1100  # at 5 °C (density of fresh water multiplied by 1.1)

# Initial position of the milk: Each of the coordinates x and y goes
# from -1 to +1 and (0, 0) is the center of the cup
x_milk, y_milk = 0, 0
# Initial depth of the milk: 0 is the surface (no milk), 1 is the bottom
z_milk = 0.4
# Initial extent of the milk: must be greater than zero
r_milk = 0.3

# Size of coffee cup (in centimeter, makes a volume of 250 mL)
Lx = 5
Ly = 5
Lz = 10

# Gravity acceleration in cm/s²
g = 980.7


# Get the default parameters, then modify them as needed
param = UserParameters()

# Choose a name for the output file
param.IO["expname"] = "coffee_with_milk"

# Select the physical quantities to save in the history file
param.IO["variables_in_history"] = ["b", "ke", "u"]

# Set the time interval (in seconds) at which these quantities are saved
param.IO["timestep_history"] = 0.01

# Make Nyles aware that CGS units are used here
param.physics["unit_length"] = "cm"
param.physics["unit_duration"] = "s"

# Turn the live animation off (or on)
param.animation["show"] = False

# Set the total length of the simulation in seconds
param.time["tend"] = 5.0
# Select automatic (or fix) time steps
param.time["auto_dt"] = True
# The following parameter is used for fix time steps
param.time["dt"] = 0.01
# The following parameters are used for automatic time steps
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.01

# Set the domain size and resolution
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz
param.discretization["global_nx"] = 32
param.discretization["global_ny"] = 32
param.discretization["global_nz"] = 64


# Initialize Nyles with these parameters
nyles = Nyles(param)

# Get access to buoyancy and its coordinates
b = nyles.model.state.b.view("i")
x = nyles.grid.x_b.view("i") / Lx
y = nyles.grid.y_b.view("i") / Ly
z = nyles.grid.z_b.view("i") / Lz

## Set the initial state of the simulation
# Define the radial coordinate r of horizontal polar coordinates around
# the point where the milk is poured into
r = np.sqrt((x-x_milk/2-1/2)**2 + (y-y_milk/2-1/2)**2)
# Choose a reference density
rho_0 = rho_coffee
# Fill the cup with coffee
b[...] = -g * rho_coffee / rho_0
# Add milk
b[z >= 1 - z_milk * np.exp(-r**2/r_milk**2)] = -g * rho_milk / rho_0

# Optionally add stirring (TODO)

# Optionally add noise to make the situation unsymmetric
# noise = np.random.normal(size=b.shape)
# b *= noise * 1e-2 + 1

# Start the simulation
nyles.run()
