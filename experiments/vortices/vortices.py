"""Experiments with quasi 2D vortices to check that Nyles works consistently."""

# Standard library imports
from enum import Enum

# Third party imports
import numpy as np

# Local Nyles imports
from nyles import Nyles


class ET(Enum):
    """Experiment Types.

    The value of each entry of this enum is a short identifier that is
    used in the experiment name.
    """
    # Rotation around the z-axis (vorticity only in k-direction)
    Vor_const = "w_const"  # constant vorticity profile: rigid body rotation
    Vor_linear = "w_lin"  # linear vorticity profile
    Vor_Gauss = "w_norm"  # Gauss-shaped vorticity profile
    Dipole = "dipole"  # two vortices of opposite sign in the horizontal plane
    TwoVortices = "2_vors"  # two vortices of the same sign
    # Rotation around the x-axis (vorticity only in i-direction)
    VerticalVor_const = "w_const_z"  # constant vorticity profile
    VerticalDipole = "dipole_z"  # two vortices of opposite sign in y,z-plane
    VerticalTwoVortices = "2_vors_z"  # two vortices of same sign


# Select a type of experiment
exp_type = ET.Vor_const

# Set size and resolution of the domain
nx = 64
ny = 64
nz = 64
Lx = 1.0
Ly = 1.0
Lz = 1.0

# Choose number of cores used
npx = 1
npy = 1
npz = 1

# Set user-defined parameters
param = {
    # I/O options
    "datadir": "~/data/Nyles",
    "expname": "vortex_{}".format(exp_type),
    "timestep_history": 1.0,
    "variables_in_history": ["vor"],
    "mode": "overwrite",
    # General model options
    "modelname": "LES",
    "geometry": "closed",
    # Grid options
    "nx": nx,
    "ny": ny,
    "nz": nz,
    "nh": 0,
    # Timestepping options
    "timestepping": "LFAM3",
    "tend": 10.0,
    "cfl": 1.5,
    # Spatial discretization options
    "orderVF": 5,  # upwind-order for vortex-force term
    "orderA": 5,  # upwind-order for advection term
}


# Initialize Nyles with the pre-defined parameters
nyles = Nyles(param)
model = nyles.model
grid = nyles.grid

# Get access to the coordinates at vorticity point
# Note that the coordinates are of Scalar type
x = grid.x_vor["i"].view("i")
y = grid.y_vor["j"].view("i")
z = grid.z_vor["k"].view("i")
# Get access to the three components of the discretized vorticity
wi = model.state.vor["i"].view("i")
wj = model.state.vor["j"].view("i")
wk = model.state.vor["k"].view("i")

### Set the initial state of the simulation
## Define the radius of the vortices
R = min([Lx, Ly]) / 2
R_yz = min([Ly, Lz]) / 2
## Define the radial coordinate r of polar coordinates (r, phi) = (x, y)
# 1) around the center (Lx/2, Ly/2)
r = np.sqrt((x-Lx/2)**2 + (y-Ly/2)**2)
# 2) shifted to the right by 1/4 of the domain length
r_right = np.sqrt((x-3*Lx/4)**2 + (y-Ly/2)**2)
# 3) shifted to the left by 1/4 of the domain length
r_left = np.sqrt((x-Lx/4)**2 + (y-Ly/2)**2)
## Define vertical polar coordinates (r, phi) = (y, z)
r_yz = np.sqrt((y-Ly/2)**2 + (z-Lz/2)**2)
r_up = np.sqrt((y-Ly/2)**2 + (z-3*Lz/4)**2)
r_down = np.sqrt((y-Ly/2)**2 + (z-Lz/4)**2)
## Define the initial vorticity field
if exp_type == ET.Vor_const:
    wk[...] = 1
    wk[r > R] = 0
elif exp_type == ET.Vor_linear:
    wk[...] = r/R
    wk[r > R] = 0
elif exp_type == ET.Vor_Gauss:
    wk[...] = np.exp(-(r/R)**2)
elif exp_type == ET.Dipole:
    wk[...] = np.exp(-(r_right/R)**2)
    wk -= np.exp(-(r_left/R)**2)
elif exp_type == ET.TwoVortices:
    wk[...] = np.exp(-(r_right/R)**2)
    wk += np.exp(-(r_left/R)**2)
elif exp_type == ET.VerticalVor_const:
    wi[...] = 1
    wi[r_yz > R_yz] = 0
elif exp_type == ET.VerticalDipole:
    wi[...] = np.exp(-(r_up/R)**2)
    wi -= np.exp(-(r_down/R)**2)
elif exp_type == ET.VerticalTwoVortices:
    wi[...] = np.exp(-(r_up/R)**2)
    wi += np.exp(-(r_down/R)**2)
else:
    raise ValueError("No experiment defined for type " + repr(exp_type))

# Optionally add noise
# for vor in [wi, wj, wk]:
#     noise = np.random.normal(size=vor.shape) * grid.msk
#     # TODO: check if these two lines are necessary
#     # noise -= grid.domain_integration(noise) * grid.msk / grid.area
#     # grid.fill_halo(noise)
#     vor += 1e-3 * noise

# Set-up the velocity field
nyles.invert_vorticity()  # TODO: implement this method

# Start the simulation
nyles.run()
