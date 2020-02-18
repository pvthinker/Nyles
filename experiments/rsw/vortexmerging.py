""" Ekman layer experiment

A uniform steady wind is imposed at the top surface of a biperiodic domain ('xy'). The Coriolis force is an O(1) term of the momentum balance. Contrary to the theory, the vertical viscosity is not prescribed. Instead, the model  resolves the full 3D turbulence sustained by the shear instability.

"""

import numpy as np

import nyles as nyles_module
import parameters
from mpi4py import MPI
import mpitools

nh = 3
factor = 4
nxglo = 32*factor
nyglo = 32*factor
nzglo = 1

npx = 1
npy = 1
npz = 1

nx = nxglo//npx
ny = nyglo//npy
nz = nzglo//npz

Lx = 1.
Ly = 1.
Lz = 0.25

# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "RSW"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "RSW_0"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['h', 'u', 'pv']

param.IO["timestep_history"] = .25  # 0.0 saves every frame
param.IO["disk_space_warning"] = 0.5  # in GB
param.IO["simplified_grid"] = True

param.time["timestepping"] = "RK3_SSP"
param.time["tend"] = 20.
param.time["auto_dt"] = False
# parameter if auto_dt is False
param.time["dt"] = 2e-3
# parameters if auto_dt is True
param.time["cfl"] = 1.2
param.time["dt_max"] = 2./factor

param.discretization["global_nx"] = nxglo
param.discretization["global_ny"] = nyglo
param.discretization["global_nz"] = nzglo
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 5  # upwind-order for advection term
param.discretization["orderKE"] = 5  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz
param.multigrid["nglue"] = 4
param.multigrid["maxite"] = 2
param.multigrid["tol"] = 1e-3


param.physics["rotating"] = True
param.physics["coriolis"] = 5.


nyles = nyles_module.Nyles(param)


h = nyles.model.state.h.view('i')
u = nyles.model.state.u['i'].view('i')
v = nyles.model.state.u['j'].view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz


def vortex(x, y, x0, y0, d):
    d2 = (x-x0)**2+(y-y0)**2
    return np.exp(-d2/(2*d**2))


def set_h(h, x_, y_):
    H = 1.
    d = 0.1
    dsep = 1.2*d
    amp = 0.05
    h[:] = vortex(x_, y_, 0.5, 0.5-dsep, d) + vortex(x_, y_, 0.5, 0.5+dsep, d)
    h *= amp
    h += H


f0 = param.physics["coriolis"]
dx = nyles.grid.dx/Lx
dy = nyles.grid.dy/Ly
set_h(h, x+0.5*dx, y+0.5*dy)
u[:,:-1,:] = -np.diff(h, axis=1)/f0*(dx/dy)
v[:,:,:-1] = +np.diff(h, axis=2)/f0*(dy/dx)
u[:,:,-1] = 0.
set_h(h, x, y)

nyles.model.diagnose_var(nyles.model.state)

nyles.run()
