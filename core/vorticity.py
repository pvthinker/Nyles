import fortran_vorticity as fortran
import variables as var
from timing import timing


@timing
def vorticity(state, fparameter):
    """
    compute the vorticity

    omega_k = delta_i[u_j] - delta_j[u_i]

    direction 'k' should be the first one,
    flipview('k') does that

    TODO:
       add the Coriolis term
       carefully handle the boundary condition: no-slip or free-slip
    """

    perm = {'i': ('k', 'j'),
            'j': ('i', 'k'),
            'k': ('j', 'i')}

    for dirk in 'ijk':
        dirj, diri = perm[dirk]

        ui = state.u[diri].flipview(dirk)
        uj = state.u[dirj].flipview(dirk)
        wk = state.vor[dirk].flipview(dirk)

        fortran.vorticity(ui, uj, wk)
        if (fparameter > 0.) and (dirk is 'k'):
            wk[:, :-1, :-1] += fparameter

        #if dirk=='j':
            # this is the future way of imposing a stress at a boundary
            # below is the implementation for the lid driven cavity
            # wk[:, :, -1] = (5e-1-uj[:,:,-1])

@timing
def vorticity2D(state, fparameter):
    """
    compute the vorticity

    omega_k = delta_i[u_j] - delta_j[u_i]

    direction 'k' should be the first one,
    flipview('k') does that

    TODO:
       carefully handle the boundary condition: no-slip or free-slip
    """

    perm = {'i': ('k', 'j'),
            'j': ('i', 'k'),
            'k': ('j', 'i')}

    for dirk in 'k':
        dirj, diri = perm[dirk]

        ui = state.u[diri].flipview(dirk)
        uj = state.u[dirj].flipview(dirk)
        wk = state.vor[dirk].flipview(dirk)

        fortran.vorticity(ui, uj, wk)
        if (fparameter > 0.) and (dirk is 'k'):
            wk[:, :-1, :-1] += fparameter

        #if dirk=='j':
            # this is the future way of imposing a stress at a boundary
            # below is the implementation for the lid driven cavity
            # wk[:, :, -1] = (5e-1-uj[:,:,-1])

def diag_rsw_pv(state):
    dirk = "k"
    vorticity = state.vor[dirk].flipview(dirk)
    h = state.h.flipview(dirk)
    pv = state.pv.flipview(dirk)
    pv[:, :-1, :-1] = h[:, :-1, :-1] + h[:, :-1, 1:] + h[:, 1:, :-1] + h[:, 1:, 1:]
    pv[:, :-1, :-1] = 4*vorticity[:, :-1, :-1]/pv[:, :-1, :-1]
            
@timing
def vorticity_all_comp(state):
    """
    purpose : reference point to compare with the previous

    verdict: it's faster!
    """

    ui = state.u['i'].view('i')
    uj = state.u['j'].view('i')
    uk = state.u['k'].view('i')

    wi = state.vor['i'].view('i')
    wj = state.vor['j'].view('i')
    wk = state.vor['k'].view('i')

    fortran.vorticity_all_comp(ui, uj, uk, wi, wj, wk)


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

    x = np.linspace(0,nx,nx)
    y = np.linspace(0,ny,ny)

    x = np.repeat(x[np.newaxis,:], ny, axis = 0)
    x = np.repeat(x[:,:,np.newaxis], nz, axis = 2)
    y = np.repeat(y[:,np.newaxis], nx, axis = 1)
    y = np.repeat(y[:,:,np.newaxis], nz, axis = 2)

    print(np.shape(u))
    print(np.shape(x))
    print(np.shape(y))

    Omega = 2

    u[...] = - Omega * y
    v[...] = Omega * x

    vorticity(state, 0.)

    vorticity_all_comp(state)

    vort = state.vor['k'].view('k')
    print(np.shape(vort))
    print(np.mean(vort)) #should be 2*omega

    plt.figure()
    plt.pcolor(vort[:,:,0])
    plt.colorbar()
    plt.show()
