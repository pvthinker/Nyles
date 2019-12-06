"""Buoyancy experiment: unstable stratification due to a temperature
gradient near the boiling point of water.

TODO: start with no stratification and build up the gradient by heating
the bottom like in a cooking pot.
"""

# Third party imports
import numpy as np
import gsw

# Local Nyles imports
from nyles import Nyles
from parameters import UserParameters


# Water temperature [°C]
T_surface = 20
T_bottom = 100
# Vertical temperature profile (linear)
T = lambda z: T_bottom + z * (T_surface - T_bottom)

# Salinity of the water [g/kg]
S = 0

# Side lenghts of a cooking pot with quadratic base area [cm]
Lx = Ly = 15
Lz = 8.89
# Volume: 15 cm * 15 cm * 8.89 cm = 2000.25 cm^3 = 2000.25 mL

# Gravity acceleration [cm/s²]
g = 981


# Get the default parameters, then modify them as needed
param = UserParameters()

# Choose a name for the output file
param.IO["expname"] = "boiling_water"

# Add two passive tracers to track the sinking of the surface layer
# and the rising of the bottom layer
param.model["n_tracers"] = 2

# Select the physical quantities to save in the history file
param.IO["variables_in_history"] = ["b", "ke", "u"]

# Set the time interval at which these quantities are saved [s]
param.IO["timestep_history"] = 0.1

# Make Nyles aware that CGS units are used here
param.physics["unit_length"] = "cm"
param.physics["unit_duration"] = "s"

# Turn the live animation for a tracer off (or use True to turn it on)
param.animation["show"] = False
param.animation["style"] = "tracer"

# Set the total length of the simulation [s]
param.time["tend"] = 10.0
# Select automatic (or fix) time steps
param.time["auto_dt"] = True
# The following parameter is used for fix time steps [s]
param.time["dt"] = 0.01
# The following parameters are used for automatic time steps
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.01

# Set the domain size
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

# Choose the spatial resolution
param.discretization["global_nx"] = 64
param.discretization["global_ny"] = 64
param.discretization["global_nz"] = 64


# Initialize Nyles with these parameters
nyles = Nyles(param)

# Get access to buoyancy, the passive tracers, and the vertical coordinate
t0 = nyles.model.state.t0.view("i")
t1 = nyles.model.state.t1.view("i")
b = nyles.model.state.b.view("i")
z = nyles.grid.z_b.view("i") / Lz

# Set the initial state of the simulation.
# We neglect the dependence of density on pressure, since this variation
# is very small in a cooking pot.
rho = gsw.rho(S, T(z), p=0)
b[...] = -g * rho / np.mean(rho)
# Add tracers to visualize the mixing
t0[z < 0.1] = 1
t1[z > 0.9] = 1

# Add noise to make the situation unsymmetric
noise = np.random.normal(size=b.shape)
b *= noise * 1e-3 + 1

# Start the simulation
nyles.run()
