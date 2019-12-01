from timing import timing

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


