"""

Projection functions to enforce div U = 0

"""


def calculate_p_from_dU(mg, state, dstate):
    """ 
    This solves the poisson
    equation with dU (dstate.U),
    stores the result in p
    (state.p) and updates dU

    mg is the multigrid object (with all data and methods)
    """

    # compute divergence
    for count, i in enumerate('jki'):
        div = state.work.view(i)
        dU = dstate.U[i].view(i)
        if count == 0: div *= 0
        div[:, :, 1:] += dU[:, :, 1:]-dU[:, :, :-1]

    # at the end of the loop div and dU are in the 'i' convention
    # this is mandatory because MG only works with the 'i' convention
        
    # copy divergence into the multigrid RHS
    # watch out, halo in MG is nh=1, it's wider for div
    b = mg.grid[0].toarray('b')
    x = mg.grid[0].toarray('x')

    # this is triplet of slices than span the MG domain (inner+MG halo)
    # typically mg_idx = (kidx, jidx, iidx)
    # with kidx = slice(k0, k1) the slice in the k direction
    # todo: implement that into state [easy once MG is back into the master branch]
    mg_idx = state.mg_idx
    b[:] = div[mg_idx]

    # solve
    mg.solve_directly()
    
    # copy MG solution to pressure
    p = state.p.view('i')
    p[mg_idx] = x[:]

    # correct dU
    # now we start with the 'i' convention
    for i in 'ijk':
        p = state.p.view(i)
        dU = dstate.U[i].view(i)
        dU[:, :, :-1] -= p[:, :, 1:]-p[:, :, :-1]
