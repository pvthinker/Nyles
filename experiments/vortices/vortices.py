"""Experiments with 2D vortices to check that Nyles works consistently."""

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
    Vor_const = "w_const"  # rigid body rotation
    Vor_linear = "w_lin"  # vortex with (limited) linear vorticity profile
    Vor_Gauss = "w_norm"  # vortex with Gauss shaped vorticity profile
    Dipole = "dipole"  # two vortices of opposite sign
    TwoVortices = "2_vors"  # two vortices of same sign
    VerticalDipole = "dipole_z"  # a dipole with rotation axes parallel to x
    VerticalTwoVortices = "2_vors_z"  # two vortices of same sign with
                                      # rotation axes parallel to x


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

# Get access to the important model variables
xr, yr, zr = grid.xr, grid.yr, grid.zr
vor = model.state.vor

### Set the initial state of the simulation
## Define the radius of the vortices
R = min([Lx, Ly]) / 2
## Define the radial coordinate r of polar coordinates (r, phi) = (x, y)
# 1) around the origin (0, 0)
r = np.sqrt(xr**2 + yr**2)
# 2) around (Lx/2, 0)
r1 = np.sqrt((xr-Lx/2)**2 + yr**2)
# 3) around (-Lx/2, 0)
r2 = np.sqrt((xr+Lx/2)**2 + yr**2)
## Define vertical polar coordinates (r, phi) = (y, z)
r3 = np.sqrt(yr**2 + (zr-Lz/2)**2)
r4 = np.sqrt(yr**2 + (zr+Lz/2)**2)
## Define the initial vorticity field
if exp_type == ET.Vor_const:
    vor[...] = 1
    vor[r > R] = 0
elif exp_type == ET.Vor_linear:
    vor[...] = r/R
    vor[r > R] = 0
elif exp_type == ET.Vor_Gauss:
    vor[...] = np.exp(-(r/R)**2)
elif exp_type == ET.Dipole:
    vor[...] = np.exp(-(r1/R)**2)
    vor -= np.exp(-(r2/R)**2)
elif exp_type == ET.TwoVortices:
    vor[...] = np.exp(-(r1/R)**2)
    vor += np.exp(-(r2/R)**2)
elif exp_type == ET.VerticalDipole:
    vor[...] = np.exp(-(r3/R)**2)
    vor -= np.exp(-(r4/R)**2)
elif exp_type == ET.VerticalTwoVortices:
    vor[...] = np.exp(-(r3/R)**2)
    vor += np.exp(-(r4/R)**2)
else:
    raise ValueError("No experiment defined for type " + repr(exp_type))

# Optionally add noise
#noise = np.random.normal(size=vor.shape) * grid.msk
# TODO: check if these two lines are necessary
##noise -= grid.domain_integration(noise) * grid.msk / grid.area
##grid.fill_halo(noise)
#vor += 1e-3 * noise

# Start the simulation
nyles.run()
