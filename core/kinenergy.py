import fortran_kinenergy as fortran
import variables as var
from timing import timing


@timing
def kinenergy(state, grid, order):
    """Compute kinetic energy from the model state.

    In z-coordinates, the kinetic energy is
      ke = 1/2 * (u**2 / dx**2 + v**2 / dy**2 + w**2 / dz**2).
    """
    for direction in 'ijk':
        # Get array of covariant velocity
        u = state.u[direction].view(direction)
        # In sigma-coordinates, use the contravariant velocity
        # U = state.U[direction].view(direction)
        # Get array of kinetic energy
        ke = state.ke.view(direction)
        if direction == "i": ke[...] = 0.
        # Get the metric term 1/dx**2 or 1/dy**2 or 1/dz**2
        ids2 = grid.ids2[direction]

        fortran.kin(u, u, ke, ids2, order)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import numpy as np

    import topology as topo
    from grid import Grid


    procs = [1, 1, 1]
    topo.topology = "closed"
    myrank = 0

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    nx = 128
    ny = 128
    nz = 128

    Lx = 1.0
    Ly = 1.0
    Lz = 1.0

    param = {
        'nx': nx, 'ny': ny, 'nz': nz,
        'Lx': Lx, 'Ly': Ly, 'Lz': Lz,
        'nh': 2, 'neighbours': neighbours,
    }
    state = var.get_state(param)
    grid = Grid(param)

    u0 = 0.8
    v0 = -4
    w0 = 0.6

    u = state.u['i'].view()
    v = state.u['j'].view()
    w = state.u['k'].view()

    u[...] = u0
    v[...] = v0
    w[...] = w0

    kinenergy(state, grid)

    ke = state.ke.view('k')

    print("max value:", np.max(ke))
    print("should be:", 1/2 * (u0**2/(Lx/nx)**2 + v0**2/(Ly/ny)**2 + w0**2/(Lz/nz)**2))

    plt.figure()
    plt.pcolor(ke[:,:,10])
    plt.colorbar()

    plt.show()
