import vortex_force as vortf
import variables as var
import tracer as tracer
import timescheme as ts
import vorticity as vort
import bernoulli as bern
import kinenergy as kinetic
import projection
import topology as topo
from timing import timing
import mg
import pickle

import numpy as np
import matplotlib.pyplot as plt

"""
LES model

At each step n, this model doesn't suppose that the velocity at step n-1 had 0 divergence.
It follows the procedure from the Ferziger p.180

"""


class LES(object):

    def __init__(self, param, grid, linear=False):
        self.nonlinear = not linear
        self.grid = grid
        self.state = var.get_state(param)
        self.timescheme = ts.Timescheme(param, self.state)
        self.timescheme.set(self.rhs, self.diagnose_var)
        self.traclist = ['b']
        self.orderA = param["orderA"]
        self.orderVF = param["orderVF"]
        self.rotating = param["rotating"]
        if self.rotating:
            # convert Coriolis parameter into its covariant quantity
            # i.e. multiply with cell horizontal area
            area = self.grid.dx*self.grid.dy
            self.fparameter = param["coriolis"] * area
        else:
            self.fparameter = 0.
        self.mg = mg.Multigrid(param)
        self.stats = []

    @timing
    def diagnose_var(self, state):
        # Diagnostic variables
        projection.calculate_p_from_dU(self.mg, state, state, self.grid)
        div = self.state.work
        projection.compute_div(div, self.state, self.grid, timing=False)
        self.update_stats()
        U_from_u(state, self.grid)
        if self.nonlinear:
            vort.vorticity(state, self.fparameter)
            kinetic.kinenergy(state, self.grid)

    @timing
    def rhs(self, state, t, dstate):
        reset_state(dstate)
        # TODO: if this function call stays here, the flag in rhstrac
        # can be removed.  Other possibility: remove the reset_state and
        # add the reset to the vortex_force term.
        # buoyancy
        tracer.rhstrac(state, dstate, self.traclist, self.orderA)

        # vortex force
        if self.nonlinear:
            vortf.vortex_force(state, dstate, self.orderVF)
        # bernoulli
        bern.bernoulli(state, dstate, self.grid)

    @timing
    def forward(self, t, dt):
        self.timescheme.forward(self.state, t, dt)
    
    def update_stats(self):
        stats = self.mg.stats
        
        div = self.state.work
        maxdiv = np.max(np.abs(div.view()))
        stats['maxdiv'] = maxdiv
        if hasattr(self, 'stats'):
            self.stats += [stats]
        else:
            self.stats = [stats]

    def write_stats(self, path):
        fid = open('%s/stats.pkl' % path, 'bw')
        pickle.dump(self.stats, fid)
        
        
@timing
def U_from_u(state, grid):
    # copied from lotsofstuff
    # for now implements only cartesian
    metric = 'cartesian'  # dx, dy and dz are uniform, though not necessarily equal

    if metric == 'cartesian':
        u = state.u['i'].view()
        v = state.u['j'].view()
        w = state.u['k'].view()

        U = state.U['i'].viewlike(state.u['i'])
        V = state.U['j'].viewlike(state.u['j'])
        W = state.U['k'].viewlike(state.u['k'])

        U[:] = u * grid.idx2
        V[:] = v * grid.idy2
        W[:] = w * grid.idz2

    elif metric == 'sigma1D':
        raise NotImplementedError
        # cf Roullet et al, OM2017
        # slope = dz/dx, at cell center
        u = state.u['i'].view('j')
        v = state.u['j'].view('j')
        w = state.u['k'].view('j')

        U = state.U['i'].view('j')
        V = state.U['j'].view('j')
        W = state.U['k'].view('j')

        V[:] = v * grid.idy2
        for j in range(ny):
            U[j][:, :] = u[j][:, :] - sxp(slope[j][:, :]*sym(w[j][:, :]))
            W[j][:, :] = gamma[j][:, :]*w[j][:, :] - \
                syp(slope[j][:, :]*sxm(u[j][:, :]))
    else:
        raise ValueError("unknown metric")


@timing
def reset_state(state):
    for var_name, var_type in state.toc.items():
        if var_type == "scalar":
            var = state.get(var_name).view()
            var *= 0.0
        else:
            for i in "ijk":
                var = state.get(var_name)[i].view()
                var *= 0.0


if __name__ == "__main__":
    from grid import Grid


    procs = [1, 1, 1]
    topo.topology = 'closed'
    myrank = 0
    nz = 32
    nh = 2

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    param = {
        'nx': 48, 'ny': 64, 'nz': nz, 'nh': nh,
        'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,
        'orderA': 5, 'orderVF': 5,
        'timestepping': 'LFAM3',
        'neighbours': neighbours,
        'procs': procs, 'topology': topo.topology,
        'npre': 3, 'npost': 3, 'omega': 0.8, 'ndeepest': 20, 'maxite': 20, 'tol': 1e-6
    }

    grid = Grid(param)

    t = 0.0
    dt = 0.1

    model = LES(param, grid)

    u = model.state.u['i'].view()
    u[:, :, :] = .1

    for kt in range(10):
        model.forward(t, dt)
        t += dt

    u = model.state.u['i'].view('k')

    plt.figure()
    plt.pcolor(u[0, :, :])
    plt.colorbar()

    plt.show()
