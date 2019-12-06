"""Buoyancy experiment: unstable stratification due to a temperature
gradient near the freezing point of water.

TODO: start with no stratification and build up the gradient by cooling
the surface like in polar regions.  Optionally add salinity flux due to
brine rejection.
"""

# Third party imports
import numpy as np
import gsw

# Local Nyles imports
from nyles import Nyles
from parameters import UserParameters


# Water temperature [°C]
T_surface = 0
T_bottom = 10
# Vertical temperature profile (linear)
T = lambda z: T_bottom + z * (T_surface - T_bottom)

# Salinity of the water [g/kg]
S_surface = 0  # increase to have a saline layer on top
S_bottom = 0
# Vertical salinity profile (thin saline layer on top)
S = lambda z: S_surface * (z > 0.95) + S_bottom * (z <= 0.95)


# Side lenghts of the basin [m]
Lx = Ly = 100
Lz = 100


# Latitude of the basin [degree]
lat = -80
# The Weddell Sea is around 80° South

# Get the default parameters, then modify them as needed
param = UserParameters()

# Choose a name for the output file
param.IO["expname"] = "freezing_water"

# Add a passive tracer to track the sinking of the surface layer
param.model["n_tracers"] = 1

# Select the physical quantities to save in the history file
param.IO["variables_in_history"] = ["b", "ke", "u"]

# Set the time interval at which these quantities are saved [s]
param.IO["timestep_history"] = 10.0

# Turn the live animation for a tracer off (or use True to turn it on)
param.animation["show"] = False
param.animation["style"] = "tracer"

# Set the total length of the simulation [s]
param.time["tend"] = 1000.0
# Select automatic (or fix) time steps
param.time["auto_dt"] = True
# The following parameter is used for fix time steps [s]
param.time["dt"] = 1.0
# The following parameters are used for automatic time steps
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.3

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

# Get access to buoyancy, the passive tracer, and the vertical coordinate
t0 = nyles.model.state.t0.view("i")
b = nyles.model.state.b.view("i")
z = nyles.grid.z_b.view("i") / Lz

# Set the initial state of the simulation.
p = gsw.p_from_z(z, lat)
g = gsw.grav(lat, p)
rho = gsw.rho(S(z), T(z), p)
b[...] = -g * rho / np.mean(rho)
# Add a tracer to visualize the mixing
t0[(z > 0.7) & (z < 0.8)] = 1

# Add noise to make the situation unsymmetric
noise = np.random.normal(size=b.shape)
b *= noise * 1e-5 + 1

# Start the simulation
nyles.run()
