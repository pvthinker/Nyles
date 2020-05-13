""" Ekman layer experiment

A uniform steady wind is imposed at the top surface of a biperiodic domain ('xy'). The Coriolis force is an O(1) term of the momentum balance. Contrary to the theory, the vertical viscosity is not prescribed. Instead, the model  resolves the full 3D turbulence sustained by the shear instability.

"""
import numpy as np

import nyles as nyles_module
import parameters
from mpi4py import MPI

nh = 3
nxglo = 32*2
nyglo = 32*2
nzglo = 32

npx = 2
npy = 2
npz = 1

nx = nxglo//npx
ny = nyglo//npy
nz = nzglo//npz

Lx = 2.
Ly = 2.
Lz = 1.

# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "perio_xy"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "ekman_2"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['b', 'u', 'vor', 'div']

param.IO["timestep_history"] = 0.5  # 0.0 saves every frame
param.IO["disk_space_warning"] = 0.5  # in GB
param.IO["simplified_grid"] = True

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 10.
param.time["auto_dt"] = False
# parameter if auto_dt is False
param.time["dt"] = 0.02
# parameters if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.8

param.discretization["global_nx"] = nxglo
param.discretization["global_ny"] = nyglo
param.discretization["global_nz"] = nzglo
param.discretization["orderVF"] = 3  # upwind-order for vortex-force term
param.discretization["orderA"] = 3  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz
param.multigrid["nglue"] = 1
param.multigrid["tol"] = 1e-4


param.physics["forced"] = True
param.physics["rotating"] = True
param.physics["coriolis"] = .25


class Forcing(object):
    def __init__(self, param, grid):
        if param["loc"][0] == param["procs"][0]-1:
            self.toplevel = True
        else:
            self.toplevel = False

        x = grid.x_u_1D / param["Lx"]
        y = grid.y_u_1D / param["Ly"]-0.5

        yy, xx = np.meshgrid(y, x, indexing="ij")

        self.tau = 1e-2#*np.exp(-yy**2/(2*.2**2))  # wind stress
        self.u_wind = 0.1
        
        self.tau *= 1./grid.dz
        self.tau *= grid.dx  # transform into covariant
        self.u_wind *= grid.dx
        self.cff = 1./grid.dz

    def add(self, state, dstate, time):
        du = dstate.u["i"].view("i")
        u = state.u["i"].view("i")
        if self.toplevel:
#            du[-1, :, :] += (self.u_wind-u[-1, :, :])*self.cff
            du[-1, :, :] += self.tau
            du[-2, :, :] += self.tau*.5


nyles = nyles_module.Nyles(param)

# the user must attach the forcing to the model
nyles.model.forcing = Forcing(nyles.param, nyles.grid)

b = nyles.model.state.b.view('i')
u = nyles.model.state.u['i'].view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz


# linear stratification
#b[:] = 1e-3*(z-0.5)

np.random.seed(nyles.myrank)
b += 1e-2*np.random.normal(0, 1, size=np.shape(b))

nyles.model.diagnose_var(nyles.model.state)

nyles.run()
