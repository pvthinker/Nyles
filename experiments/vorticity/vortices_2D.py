"""Experiments with quasi 2D vortices to check that Nyles works consistently."""

# Standard library imports
from enum import Enum

# Third party imports
import numpy as np

# Local Nyles imports
from nyles import Nyles
from parameters import UserParameters


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


# Get the default parameters, then modify them as needed
param = UserParameters()

# Choose a name for the experiment
param.IO["expname"] = "vortex_{}".format(exp_type)

# Set the size of the domain
param.model["Lx"] = 1.0
param.model["Ly"] = 1.0
param.model["Lz"] = 1.0

# Choose the resolution of the domain
param.discretization["global_nx"] = 32
param.discretization["global_ny"] = 32
param.discretization["global_nz"] = 32

# Set the total length (in time) of the simulation
param.time["tend"] = 20.0

# Choose number of cores used
param.MPI["npx"] = 1
param.MPI["npy"] = 1
param.MPI["npz"] = 1

# Set halo size
param.MPI["nh"] = 3


# Initialize Nyles
nyles = Nyles(param)
model = nyles.model
grid = nyles.grid

# Get access to the three components of vorticity and their (x,y,z)-coordinates
wi = model.state.vor["i"].view("i")
x_wi = grid.x_vor["i"].view("i")
y_wi = grid.y_vor["i"].view("i")
z_wi = grid.z_vor["i"].view("i")

wj = model.state.vor["j"].view("i")
x_wj = grid.x_vor["j"].view("i")
y_wj = grid.y_vor["j"].view("i")
z_wj = grid.z_vor["j"].view("i")

wk = model.state.vor["k"].view("i")
x_wk = grid.x_vor["k"].view("i")
y_wk = grid.y_vor["k"].view("i")
z_wk = grid.z_vor["k"].view("i")

### Set the initial state of the simulation
## Define the radius of the vortices
R_wk = min([grid.Lx, grid.Ly]) / 4
R_wi = min([grid.Ly, grid.Lz]) / 4
## Define the radial coordinate r of polar coordinates (r, phi) = (x, y)
# 1) around the center (Lx/2, Ly/2)
r_wk = np.sqrt((x_wk-grid.Lx/2)**2 + (y_wk-grid.Ly/2)**2)
# 2) shifted to the right by 1/4 of the domain length
r_right_wk = np.sqrt((x_wk-3*grid.Lx/4)**2 + (y_wk-grid.Ly/2)**2)
# 3) shifted to the left by 1/4 of the domain length
r_left_wk = np.sqrt((x_wk-grid.Lx/4)**2 + (y_wk-grid.Ly/2)**2)
## Define vertical polar coordinates (r, phi) = (y, z)
r_wi = np.sqrt((y_wi-grid.Ly/2)**2 + (z_wi-grid.Lz/2)**2)
r_up_wi = np.sqrt((y_wi-grid.Ly/2)**2 + (z_wi-3*grid.Lz/4)**2)
r_down_wi = np.sqrt((y_wi-grid.Ly/2)**2 + (z_wi-grid.Lz/4)**2)
## Define the initial vorticity field
if exp_type == ET.Vor_const:
    wk[...] = 1
    wk[r_wk > R_wk] = 0
elif exp_type == ET.Vor_linear:
    wk[...] = r_wk/R_wk
    wk[r_wk > R_wk] = 0
elif exp_type == ET.Vor_Gauss:
    wk[...] = np.exp(-(r_wk/R_wk)**2)
elif exp_type == ET.Dipole:
    wk[...] = np.exp(-(r_right_wk/R_wk)**2)
    wk -= np.exp(-(r_left_wk/R_wk)**2)
elif exp_type == ET.TwoVortices:
    wk[...] = np.exp(-(r_right_wk/R_wk)**2)
    wk += np.exp(-(r_left_wk/R_wk)**2)
elif exp_type == ET.VerticalVor_const:
    wi[...] = 1
    wi[r_wi > R_wi] = 0
elif exp_type == ET.VerticalDipole:
    wi[...] = np.exp(-(r_up_wi/R_wi)**2)
    wi -= np.exp(-(r_down_wi/R_wi)**2)
elif exp_type == ET.VerticalTwoVortices:
    wi[...] = np.exp(-(r_up_wi/R_wi)**2)
    wi += np.exp(-(r_down_wi/R_wi)**2)
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
