import vortex_force as vortf
import variables as var
import tracer
import timescheme as ts
import vorticity as vort
import bernoulli as bern
import kinenergy as kinetic
import viscosity as visc
import projection
import boundarycond as bc
import topology as topo
from timing import timing
import mg
import pickle
import cov_to_contra
import halo

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
        self.traclist = ['b']
        # Add tracers (if any)
        for i in range(param["n_tracers"]):
            t_nickname = "t{}".format(i)
            t_name = "tracer{}".format(i)
            self.traclist.append(t_nickname)
            var.modelvar[t_nickname] = var.ModelVariable(
                'scalar', t_name,  dimension='', prognostic=True)
        self.state = var.get_state(param)
        self.halo = halo.set_halo(param, self.state)
        self.neighbours = param["neighbours"]
        self.timescheme = ts.Timescheme(param, self.state)
        self.timescheme.set(self.rhs, self.diagnose_var)
        self.orderA = param["orderA"]
        self.orderVF = param["orderVF"]
        self.rotating = param["rotating"]
        self.forced = param["forced"]
        self.diff_coef = param['diff_coef']
        self.add_viscosity = "u" in self.diff_coef.keys()
        if self.add_viscosity:
            self.viscosity = self.diff_coef['u']

        self.tracer = tracer.Tracer_numerics(
            grid, self.traclist, self.orderA, self.diff_coef)
        if self.rotating:
            # convert Coriolis parameter (in s^-1) into its covariant quantity
            # i.e. multiply with cell horizontal area
            area = self.grid.dx*self.grid.dy
            self.fparameter = param["coriolis"] * area
        else:
            self.fparameter = 0.
        self.mg = mg.Multigrid(param, grid)
        self.stats = []

    @timing
    def diagnose_var(self, state):
        #bc.apply_bc_on_velocity(state, self.neighbours)
        self.halo.fill(state.b)
        self.halo.fill(state.u)
        # Diagnostic variables
        cov_to_contra.U_from_u(state, self.grid)
        projection.compute_p(self.mg, state, self.grid)
        self.halo.fill(state.u)
        cov_to_contra.U_from_u(state, self.grid)

        # this computation is only to check the divergence
        # after the projection, this could be drop if
        # we don't need to know the information (that can
        # always be estimated offline, from 'u')
        projection.compute_div(self.state, timing=False)
        self.halo.fill(state.div)
        self.update_stats()

        if self.nonlinear:
            vort.vorticity(state, self.fparameter)
            #bc.apply_bc_on_vorticity(state, self.neighbours)
            kinetic.kinenergy(state, self.grid)
            self.halo.fill(state.vor)
            self.halo.fill(state.ke)

    @timing
    def rhs(self, state, t, dstate, last=False):
        reset_state(dstate)
        # TODO: if this function call stays here, the flag in rhstrac
        # can be removed.  Other possibility: remove the reset_state and
        # add the reset to the vortex_force term.
        # buoyancy
        self.tracer.rhstrac(state, dstate)
        # vortex force
        if self.nonlinear:
            vortf.vortex_force(state, dstate, self.orderVF)
        # bernoulli
        bern.bernoulli(state, dstate, self.grid)

        if last and self.add_viscosity:
            visc.add_viscosity(self.grid, state, dstate, self.viscosity)

        if self.forced:
            self.forcing.add(state, dstate, t)

    @timing
    def forward(self, t, dt):
        self.timescheme.forward(self.state, t, dt)
        return self.mg.stats['blowup']

    def update_stats(self):
        stats = self.mg.stats

        div = self.state.div
        maxdiv = np.max(np.abs(div.view()))
        stats['maxdiv'] = maxdiv
        if hasattr(self, 'stats'):
            self.stats += [stats]
        else:
            self.stats = [stats]

    def write_stats(self, path):
        fid = open('%s/stats.pkl' % path, 'bw')
        pickle.dump(self.stats, fid)


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
        'rotating': True, 'coriolis': 1.0,
        'neighbours': neighbours,
        'procs': procs, 'myrank': myrank,
        'npre': 3, 'npost': 3, 'omega': 0.8, 'ndeepest': 20, 'maxite': 20, 'tol': 1e-6,
        'nglue': 16, 'ncellscoarsest': 32,
        'n_tracers': 0,
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
