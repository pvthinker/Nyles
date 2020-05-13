import numpy as np

import nyles as nyles_module
import parameters


nh = 3
nxglo = 32
nyglo = 32
nzglo = 16

npx = 1
npy = 1
npz = 1

nx = nxglo//npx
ny = nyglo//npy
nz = nzglo//npz

Lx = 4e3
Ly = 4e3
Lz = 2e3

# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "forced_plume_nudging"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['b', 'u', 'vor', 'div']

param.IO["timestep_history"] = 600.  # 0.0 saves every frame
param.IO["disk_space_warning"] = 0.5  # in GB
param.IO["simplified_grid"] = True

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 3600.*24
param.time["auto_dt"] = True
# parameter if auto_dt is False
param.time["dt"] = 200.
# parameters if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = param.time["dt"]

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
param.physics["coriolis"] = 1e-4

def stratif(z):
    return 1e-2*(z-0.5)


class Forcing(object):
    def __init__(self, param, grid):
        x = grid.x_b.view('i') / param["Lx"]-0.5
        y = grid.y_b.view('i') / param["Ly"]-0.5
        z = grid.z_b.view('i') / param["Lz"]
        d = np.sqrt( x**2+y**2)
        msk = 0.5*(1.-np.tanh(d/0.1))
        self.Q = 1e-5*np.exp(-z/0.02)*msk
        self.bclim = stratif(z)
        d0 = 0.4
        width = 0.05
        dampingcoef = 1./200
        self.damping = 0.5*(1+np.tanh( (d-d0)/width ))
        self.damping *= dampingcoef
        
    def add(self, state, dstate, time):
        db = dstate.b.view("i")
        b = state.b.view("i")
        db += self.Q - self.damping*(b-self.bclim)

nyles = nyles_module.Nyles(param)

# the user must attach the forcing to the model
nyles.model.forcing = Forcing(nyles.param, nyles.grid)

b = nyles.model.state.b.view('i')
u = nyles.model.state.u['i'].view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz

# linear stratification
b[:] = stratif(z)

nyles.model.diagnose_var(nyles.model.state)

nyles.run()
