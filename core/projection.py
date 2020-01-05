"""

Projection functions to enforce div U = 0

"""
import numpy as np
from timing import timing
import boundarycond as bc

@timing
def compute_div(state, **kwargs):
    """Compute divergence from U, the contravariant velocity

    div = delta[U*vol]/vol

    here we assume that vol is uniform (Cartesian coordinates case)
    so

    div = delta[U]
    """
    for count, i in enumerate('jki'):
        div = state.div.view(i)
        dU = state.U[i].view(i)
        if count == 0:
            div[:, :, 0] = dU[:, :, 0]
            div[:, :, 1:] = np.diff(dU)
        else:
            div[:, :, 0] += dU[:, :, 0]
            div[:, :, 1:] += np.diff(dU)


@timing
def compute_p(mg, state, grid, ngbs):
    """
    This solves the poisson equation with U (state.U),
    stores the result in p (state.p)
    and corrects u (state.u) with -delta[p]

    note that the real pressure is p/dt

    mg is the multigrid object (with all data and methods)

    grid is the Grid object with the metric tensor
    """
    div = state.div
    compute_div(state)

    # at the end of the loop div and U are in the 'i' convention
    # this is mandatory because MG only works with the 'i' convention

    # copy divergence into the multigrid RHS
    # watch out, halo in MG is nh=1, it's wider for div
    b = mg.grid[0].toarray('b')
    x = mg.grid[0].toarray('x')

    # this is triplet of slices than span the MG domain (inner+MG halo)
    # typically mg_idx = (kidx, jidx, iidx)
    # with kidx = slice(k0, k1) the slice in the k direction
    mg_idx = state.div.mg_idx
    b[:] = div.view('i')[mg_idx]
#    x[:] = 0.
    # solve
    mg.solve_directly()

    # copy MG solution to pressure
    p = state.p.view('i')
    p[mg_idx] = x[:]
    #bc.apply_bc_on_p(state, ngbs)

    # correct u (the covariant component)
    # now we start with the 'i' convention
    for i in 'ijk':
        p = state.p.view(i)
        u = state.u[i].view(i)
        u[:, :, :-1] -= np.diff(p)


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    import topology as topo
    import mg
    import variables as var
    import vorticity as vort
    import grid as grid_module

    procs = [1, 1, 1]
    topo.topology = 'closed'
    myrank = 0
    nx, ny, nz = 48, 24, 32
    nh = 2

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    param = {
        'nx': nx, 'ny': ny, 'nz': nz, 'nh': nh,
        'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,
        'timestepping': 'LFAM3',
        'neighbours': neighbours,
        'procs': procs, 'myrank': myrank,
        'npre': 3, 'npost': 3, 'omega': 0.8, 'ndeepest': 20, 'maxite': 20, 'tol': 1e-12
    }

    grid = grid_module.Grid(param)

    s = var.get_state(param)
    ds = s.duplicate_prognostic_variables()
    mg = mg.Multigrid(param)
    u = ds.u['i'].view('i')
    v = ds.u['j'].view('i')
    w = ds.u['k'].view('i')

    print('set up a localized jet along the i-direction')
    k1 = nz//3
    k2 = 2*nz//3
    j1 = ny//3
    j2 = 2*ny//3
    u[k1:k2, j1:j2, :-1] = 1.
    v[:, :-1, :] = 0
    w[:-1, :, :] = 0

    U = s.U['i'].view('i')
    V = s.U['j'].view('i')
    V = s.U['k'].view('i')
    U[:] = u*grid.idx2
    V[:] = v*grid.idy2
    W[:] = w*grid.idz2

    print('and make it divergent-free')
    compute_p(mg, s, ds, grid)

    p = s.p.view('i')

    div = s.div
    compute_div(ds)

    u = ds.u['i'].view('i')
    v = ds.u['j'].view('i')
    w = ds.u['k'].view('i')

    plt.close('all')
    plt.ion()

    d = div.view('i')
    fields = [('div', d), ('u', u), ('v', v), ('w', w)]
    fig, ax = plt.subplots(2, 2)
    fig.set_size_inches([15,  8])
    for i in range(4):
        a = ax[i//2, i % 2]
        name, f = fields[i]
        im = a.imshow(f[10, :, :], origin='xy')
        a.set_xlabel('x')
        a.set_ylabel('y')
        a.set_title(name)
        plt.colorbar(im, ax=a)

    plt.show()
    for direc in 'ijk':
        s.u[direc].view('i')[:] = ds.u[direc].view('i')[:]

    vort.vorticity(s, 0.)
    vor = s.vor
    wy = vor['j'].view('i')
    plt.figure()
    plt.imshow(wy[:, :, 10], origin='xy')
    plt.colorbar()
