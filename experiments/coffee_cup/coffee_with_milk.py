"""Buoyancy experiment: mix hot coffee with cold milk."""

# Third party imports
import numpy as np

# Local Nyles imports
from nyles import Nyles


# Densities in kg m^-3
rho_coffee = 970  # at 85 °C (density of fresh water)
rho_milk = 1100  # at 5 °C (density of fresh water multiplied by 1.1)

# Initial position of the milk: Each of the coordinates x and y goes
# from -1 to +1 and (0, 0) is the center of the cup
x_milk, y_milk = 0, 0
# Initial depth of the milk: 0 is the surface (no milk), 1 is the bottom
z_milk = 0.5
# Initial extent of the milk: must be greater than zero
r_milk = 0.2

# Size of coffee cup (in centimeter, makes a volume of 225 ml)
Lx = 5.0
Ly = 5.0
Lz = 10.0

# Gravity acceleration in m/s²
g = 9.81

# Choose resolution
nx = 8
ny = 8
nz = 16
# Choose number of CPU cores used
npx = 1
npy = 1
npz = 1

# Set user-defined parameters
param = {
    # I/O options
    "datadir": "~/data/Nyles",
    "expname": "coffee_with_milk",
    "timestep_history": 1.0,
    "variables_in_history": ["b", "ke", "u"],
    "mode": "count",
    # General model options
    "modelname": "LES",
    "geometry": "closed",
    # Domain options
    "Lx": Lx,
    "Ly": Ly,
    "Lz": Lz,
    # Grid options
    "nx": nx,
    "ny": ny,
    "nz": nz,
    "nh": 0,
    # Timestepping options
    "timestepping": "LFAM3",
    "tend": 20.0,
    "cfl": 1.5,
    # Spatial discretization options
    "orderVF": 5,  # upwind-order for vortex-force term
    "orderA": 5,  # upwind-order for advection term
}


# Initialize Nyles with the pre-defined parameters
nyles = Nyles(param)
model = nyles.model
grid = nyles.grid

# Get access to buoyancy and its coordinates
b = model.state.b.view("i")
x = grid.x_b.view("i")
y = grid.y_b.view("i")
z = grid.z_b.view("i")

## Set the initial state of the simulation
# Define the radial coordinate r of horizontal polar coordinates around
# the point where the milk is poured into
r = np.sqrt((x-x_milk*Lx/2-Lx/2)**2 + (y-y_milk*Ly/2-Ly/2)**2)
# Choose a reference density
rho_0 = rho_coffee
# Fill the cup with coffee
b[...] = -g * rho_coffee/rho_0
# Add milk
b[z > Lz - Lz * z_milk * np.exp(r**2/r_milk**2)] = -g * rho_milk/rho_0

# TODO: optionally add stirring

# Optionally add noise
# noise = np.random.normal(size=b.shape) * grid.msk
# # TODO: check if these two lines are necessary
# # noise -= grid.domain_integration(noise) * grid.msk / grid.area
# # grid.fill_halo(noise)
# b += 1e-3 * noise

# Start the simulation
nyles.run()
