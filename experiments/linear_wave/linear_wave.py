"""Simulate an interface wave with linear dynamics in Nyles."""

import numpy as np

from nyles import Nyles
from parameters import UserParameters


# Set the resolution of the domain
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

param.model["modelname"] = "linear"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["expname"] = "linear_wave"
param.IO["variables_in_history"] = ["b", "u"]

# Set the length of the simulation
param.time["tend"] = 30.0
param.time["auto_dt"] = True
# The following parameter is used if auto_dt is False
param.time["dt"] = 0.5
# The following parameters are used if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.5

param.discretization["global_nx"] = nx
param.discretization["global_ny"] = ny
param.discretization["global_nz"] = nz


nyles = Nyles(param)

# Get access to the buoyancy and its coordinates
b = nyles.model.state.b.view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz

# Define the initial wave on the interface of two layers of different density.
# The amplitude of the wave must be small compared to 1 to justify the
# assumptions of linear dynamics.
amplitude = 0.03
b[:] = 1 + np.tanh((amplitude * np.cos(4*np.pi*(x+y)) + z - 1/2) / 0.02)

# Start the simulation
nyles.run()
