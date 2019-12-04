"""Simulate an interface wave with linear dynamics in Nyles."""

# Third party imports
import numpy as np

# Local Nyles imports
from nyles import Nyles
from parameters import UserParameters


# Choose the resolution of the domain
nx = 64
ny = 64
nz = 64

# Choose the size of the domain
# It must be possible to set the following lengths to arbitrary values,
# but currently, due to a problem in the handling of the metric in the
# calculation of p, it is necessary to ensure dx = dy = dz = 1. #TODO
Lx = 1.0 * nx
Ly = 1.0 * ny
Lz = 1.0 * nz


# Get the default parameters, then modify them as needed
param = UserParameters()

# Choose the linear model
param.model["modelname"] = "linear"

# Choose a name for the output file
param.IO["expname"] = "linear_wave"

# Select the physical quantities to save in the history file
param.IO["variables_in_history"] = ["b", "u"]

# Turn the live animation on (or off)
param.animation["show"] = True
# Tell the animation module that the wave is stable
param.animation["stable_stratification"] = True
# Let the animation update after every second iteration step
param.animation["iterations_per_frame"] = 2

# Set the total length of the simulation
param.time["tend"] = 30.0
# Select automatic (or fix) time steps
param.time["auto_dt"] = True
# The following parameter is used for fix time steps
param.time["dt"] = 0.5
# The following parameters are used for automatic time steps
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.5

# Set the domain size and resolution
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz
param.discretization["global_nx"] = nx
param.discretization["global_ny"] = ny
param.discretization["global_nz"] = nz


# Initialize Nyles with these parameters
nyles = Nyles(param)

# Get access to buoyancy and its coordinates
b = nyles.model.state.b.view("i")
x = nyles.grid.x_b.view("i") / Lx
y = nyles.grid.y_b.view("i") / Ly
z = nyles.grid.z_b.view("i") / Lz

# Define the initial wave on the interface of two layers of different density.
# The amplitude of the wave must be small compared to 1 to justify the
# assumptions of linear dynamics.
amplitude = 0.03
b[:] = 1 + np.tanh((amplitude * np.cos(4 * np.pi * (x + y)) + z - 1 / 2) / 0.02)

# Start the simulation
nyles.run()
