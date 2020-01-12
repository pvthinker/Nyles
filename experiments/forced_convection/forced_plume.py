import numpy as np

import nyles as nyles_module
import parameters


nh = 3
nxglo = 64
nyglo = 32
nzglo = 64

npx = 1
npy = 1
npz = 1

nx = nxglo//npx
ny = nyglo//npy
nz = nzglo//npz

Lx = 8.
Ly = 4.
Lz = 8.

# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "forced_plume_0"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['b', 'u', 'vor', 'div']

param.IO["timestep_history"] = 1.  # 0.0 saves every frame
param.IO["disk_space_warning"] = 0.5  # in GB

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 200.0
param.time["auto_dt"] = True
# parameter if auto_dt is False
param.time["dt"] = 0.2
# parameters if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.8

param.discretization["global_nx"] = nxglo
param.discretization["global_ny"] = nyglo
param.discretization["global_nz"] = nzglo
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 5  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz
param.multigrid["nglue"] = 1


param.physics["forced"] = True
param.physics["rotating"] = True
param.physics["coriolis"] = 1.

class Forcing(object):
    def __init__(self, param, grid):
        x = grid.x_b.view('i') / param["Lx"]-0.5
        y = grid.y_b.view('i') / param["Ly"]-0.5
        z = grid.z_b.view('i') / param["Lz"]
        d = np.sqrt( x**2+y**2)
        msk = 0.5*(1.-np.tanh(d/0.1))
        self.Q = 1e-1*np.exp(-z/0.02)*msk
        
    def add(self, state, dstate, time):
        db = dstate.b.view("i")
        db += self.Q

nyles = nyles_module.Nyles(param)

# the user must attach the forcing to the model
nyles.model.forcing = Forcing(nyles.param, nyles.grid)

b = nyles.model.state.b.view('i')
u = nyles.model.state.u['i'].view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz

# linear stratification
b[:] = 1e-1*(z-0.5)

nyles.model.diagnose_var(nyles.model.state)

nyles.run()
