import fortran_kinenergy as fortran
import variables as var
from timing import timing


@timing
def kinenergy(state):
    """
    compute the kinetic energy function from the model state

    """
    # loop over direction
    for direction in 'ijk':
        # Get array of covariant velocity u
        u = state.u[direction].view(direction)
        # In sigma-coordinates, use the contravariant velocity U
        # U = state.U[direction].view(direction)

        # Get array of kinetic energy ke
        ke = state.ke.view(direction)
        # In z-coordinates, the kinetic energy is ke = u**2 / ds**2,
        # where ds2 is the diagonal term of the inverse metric tensor
        ds2 = 1.  # TODO: use dx**-2, dy**-2, dz**-2

        if direction == 'i':
            fortran.kin(u, u, ke, ds2, 1)  # overwrite ke
        else:
            fortran.kin(u, u, ke, ds2, 0)  # add to ke


# ----------------------------------------------------------------------
if __name__ == '__main__':

    import topology as topo
    import matplotlib.pyplot as plt
    import numpy as np

    procs = [1, 1, 1]
    topo.topology = "closed"
    myrank = 0

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    nx = 128
    ny = 128
    nz = 128

    param = {'nx': nx, 'ny': ny, 'nz': nz, 'nh': 2, 'neighbours' : neighbours}
    state = var.get_state(param)

    u = state.u['i'].view('k')
    v = state.u['j'].view('k')
    w = state.u['k'].view('k')

    u[...] = 4
    v[...] = 8
    w[...] = 6

    kinenergy(state)

    ke = state.ke.view('k')

    print(np.max(ke))
    print("should be : ", np.max(u)**2/2 + np.max(v)**2/2 + np.max(w)**2/2)

    plt.figure()
    plt.pcolor(ke[:,:,10])
    plt.colorbar()

    plt.show()
