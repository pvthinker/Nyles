from nyles import Nyles
from param import Param
from grid import Grid

import numpy as np

param = Param()
param.modelname = 'boussinesq'
param.expname = 'taylorinst_00'

# domain and resolution
param.nx = 64
param.ny = 64
param.nz = 64
param.npx = 1
param.npy = 1
param.npz = 1

param.Lx = 1.0
param.Ly = 1.0
param.Lz = 1.0
param.geometry = 'closed'

# time
param.tend = 10.0
param.cfl = 1.5
param.adaptable_dt = True
param.dt = 0.004
param.dtmax = 1e-2

# discretization
param.order = 5

# output
param.plot_var = 'b'
param.var_to_save = ['u', 'b', 'p']
param.list_diag = 'all'
param.freq_his = 1.0
param.freq_diag = 1.0

# plot
param.plot_interactive = True
param.freq_plot = 10
param.colorscheme = 'imposed'
param.cax = [-0.6, 0.6]
param.generate_mp4 = True

# physics
param.forcing = False
param.diffusion = False
param.noslip = False

grid = Grid(param)

xr, yr, zr = grid.xr, grid.yr, grid.zr
Lx, Ly, Lz = param.Lx, param.Ly, param.Lz

nyles = Nyles(param, grid)
model = nyles.model

xr, yr = grid.xr, grid.yr
buoy = model.var.get('b')


def sigmoid(x, delta):
    return 1/(1+np.exp(-(x-0.5)/delta))


def stratif():
    sigma = grid.dx/2  # width of the interface
    b = sigmoid(zr/param.Lz, sigma/param.Lz)
    return b


buoy[:] = (1-stratif() - 0.5)
# add noise to trigger the instability
noise = np.random.normal(size=np.shape(yr), scale=1.) * grid.msk
noise -= grid.domain_integration(noise) * grid.msk / grid.area
grid.fill_halo(noise)

buoy += 1e-3 * noise

nyles.loop()
