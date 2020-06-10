""" Rayleigh Benard experiment

A closed box is warm uniformally from the bottom with a heat flux +Q and cool uniformally from the top with a heat flux -Q.

The heat transport from bottom to top is done by convections cells.

Current limitation (feb 2020): the boundary conditions (BC) for the momentum are free-slip. It would be better to have no-slip BC. This is an on-going work in the code development.

"""
import numpy as np

import nyles as nyles_module
import parameters

Kdiff = 1e-2
visc = 1e-2

# imposed flux or imposed temperature
forcing_condition = "temp" # "flux" or "temp"

Q = 5e-4 # value of the flux in the "flux" case
Deltab = 10e-2 # delta buoyancy in the "temp" case

nh = 3
nxglo = 32
nyglo = 4
nzglo = 4

npx = 1
npy = 1
npz = 1

nx = nxglo//npx
ny = nyglo//npy
nz = nzglo//npz

# this trick imposes that all experiments have Lz=1
# and dx==dy==dz
Lz = 1.
Lx = Lz/nzglo*nxglo
Ly = Lz/nzglo*nyglo


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "closed"#perio_xy"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "RB_00"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['b', 'u']
param.IO["simplified_grid"] = True

param.IO["timestep_history"] = 100.  # 0.0 saves every frame
param.IO["disk_space_warning"] = 0.5  # in GB

# this is the maximum time step imposed by diffusion
# when the flow becomes turbulent, the time step
# might need to be reduced. This can be done by
# using auto_dt_dt = True
dz = Lz/nzglo
maxdt = 0.1*dz**2/Kdiff

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 3000.0
param.time["auto_dt"] = True
# parameter if auto_dt is False
param.time["dt"] = maxdt
# parameters if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = maxdt

param.discretization["global_nx"] = nxglo
param.discretization["global_ny"] = nyglo
param.discretization["global_nz"] = nzglo
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 5  # upwind-order for advection term
param.discretization["orderKE"] = 5 # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz
param.multigrid["nglue"] = 4
param.multigrid["maxite"] = 4
param.multigrid["tol"] = 1e-4


param.physics["forced"] = True
if param.time["timestepping"] == "RK3_SSP":
    # the 1.5=3/2 coefficient is to compensate for the weights in
    # RK3_SSP: diffusion and viscosity are added only on the
    # third stage, at n+1/2, that enters the tendency with a
    # 2/3 coefficient
    param.physics["diff_coef"] = {"b": Kdiff*1.5, "u": visc*1.5}

else:
    param.physics["diff_coef"] = {"b": Kdiff, "u": visc}


if forcing_condition == "temp":
    param.user["param1"] = Deltab

elif forcing_condition == "flux":
    param.user["param1"] = Q

else:
    raise ValueError

class Forcing_flux(object):
    def __init__(self, param, grid):
        self.toplevel = False
        self.botlevel = False
        
        if param["loc"][0] == param["npz"]-1:
            self.toplevel = True
            
        if param["loc"][0] == 0:
            self.botlevel = True

        x = grid.x_b.view('i') / param["Lx"]-0.5
        y = grid.y_b.view('i') / param["Ly"]-0.5
        z = grid.z_b.view('i') / param["Lz"]

        self.Q = param["param1"]
        self.dz = grid.dz
        
    def add(self, state, dstate, time):
        """ Q has the dimensions of [b]*(L*T^-1) """
        db = dstate.b.view("i")        
        if self.botlevel: db[0, :, :] += self.Q/self.dz
        if self.toplevel: db[-1, :, :] -= self.Q/self.dz

class Forcing_temp(object):
    def __init__(self, param, grid):
        self.toplevel = False
        self.botlevel = False
        
        if param["loc"][0] == param["npz"]-1:
            self.toplevel = True
            
        if param["loc"][0] == 0:
            self.botlevel = True

        x = grid.x_b.view('i') / param["Lx"]-0.5
        y = grid.y_b.view('i') / param["Ly"]-0.5
        z = grid.z_b.view('i') / param["Lz"]

        self.bsurf = param["param1"]/2
        self.dz = grid.dz
        Kdiff = param["diff_coef"]["b"]
        self.coef = 2*Kdiff/dz**2
        
    def add(self, state, dstate, time):
        """ Q has the dimensions of [b]*(L*T^-1) """
        db = dstate.b.view("i")
        b = state.b.view("i")
        if self.botlevel: db[0, :, :] += self.coef*(self.bsurf-b[0])
        if self.toplevel: db[-1, :, :] += self.coef*(-self.bsurf-b[-1])
        
nyles = nyles_module.Nyles(param)

# the user must attach the forcing to the model
if forcing_condition == "temp":
    nyles.model.forcing = Forcing_temp(nyles.param, nyles.grid)
elif forcing_condition == "flux":
    nyles.model.forcing = Forcing_flux(nyles.param, nyles.grid)


b = nyles.model.state.b.view('i')
u = nyles.model.state.u['i'].view('i')
x = nyles.grid.x_b.view('i')
y = nyles.grid.y_b.view('i')
z = nyles.grid.z_b.view('i')

if forcing_condition == "temp":
    b[:] = (Lz*.5-z)/Lz*Deltab

elif forcing_condition == "flux":
    b[:] = Q*(Lz*.5-z)/Kdiff

nyles.model.diagnose_var(nyles.model.state)

nyles.run()
