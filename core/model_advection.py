"""All the tools needed to advect a tracer with an imposed velocity field."""

import variables as var
import tracer
import cov_to_contra
import projection
import timescheme as ts
from timing import timing
import pickle
import halo


class Advection(object):
    """Pure advection model.

    It advects the buoyancy with a prescribed contravariant velocity.
    """

    def __init__(self, param, grid):
        self.state = var.get_state(param)
        self.traclist = ['b']
        self.order = param['orderA']
        self.halo = halo.set_halo(param, self.state)
        self.timescheme = ts.Timescheme(param, self.state)
        self.timescheme.set(self.rhs, self.diagnose_var)
        self.param = param
        self.grid = grid
        self.stats = []

    def diagnose_var(self, state):
        # no diagnostic variables in this model
        self.halo.fill(state.b)

    def make_u_divergentfree(self):
        import mg
        mg = mg.Multigrid(self.param, self.grid)

        self.halo.fill(self.state.u)

        cov_to_contra.U_from_u(self.state, self.grid)
        projection.compute_p(mg, self.state, self.grid)
        cov_to_contra.U_from_u(self.state, self.grid)

    @timing
    def forward(self, t, dt):
        self.timescheme.forward(self.state, t, dt)
        return False

    @timing
    def rhs(self, state, t, dstate):
        tracer.rhstrac(state, dstate, self.traclist, self.order)

    def update_stats(self):
        b = self.state.b
        bmean = np.mean(b.view())
        bstd = np.bstd(b.view())
        stats = {'bmean': bmean, 'bstd': bstd}

    def write_stats(self, path):
        fid = open('%s/stats.pkl' % path, 'bw')
        pickle.dump(self.stats, fid)


if __name__ == '__main__':
    import numpy as np
    from matplotlib import pyplot as plt
    import topology as topo

    procs = [1, 1, 1]
    topo.topology = 'closed'
    myrank = 0
    nh = 2

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    nz, ny, nx = 32, 32, 128
    Lx, Ly, Lz = 1.0, 1.0, 1.0
    dx, dy, dz = Lx/nx, Ly/ny, Lz/nz
    param = {
        'nx': nx, 'ny': ny, 'nz': nz, 'nh': nh,
        'Lx': Lx, 'Ly': Ly, 'Lz': Lz,
        'neighbours': neighbours,
        # Choose a timestepping method
        'timestepping': 'LFAM3',
        # 'timestepping': 'EF',
        # Set the order of the upwind scheme
        'orderA': 5,
        # 'orderA': 3,
        # 'orderA': 1,
    }

    model = Advection(param)

    # Set up a uniform velocity along i.
    # This does not enforce the no-flow BC,
    # but that's not a problem if we don't
    # integrate too long.
    U0 = 1.0
    Ui = model.state.U['i'].view()
    Ui[...] = U0

    # Set up a Gaussian shape of the buoyancy along x
    b = model.state.b.view('i')
    for i in range(nx):
        x = (i+0.5)*dx-0.2*Lx
        b[:, :, i] = np.exp(-0.5*(x/(dx*5))**2)

    # Integrate the model forward in time and show the result
    plt.figure()
    t = 0.0
    cfl = 0.5
    dt = cfl * U0
    for kt in range(41):
        model.forward(t, dt)
        t += dt
        if kt % 5 == 0:
            plt.plot(b[0, 0, :], label='%i' % kt)
    plt.legend()
    plt.show()
