import numpy as np

def apply_bc_on_velocity(state, ngbs):
    """ apply no-flow BC on left walls, if any """
    if (-1, 0, 0) in ngbs.keys():
        pass
    else:
        w = state.u['k'].view('k')
        w[:, :, 0] = 0.
        v = state.u['j'].view('k')
        u = state.u['i'].view('k')
        v[:,:,0] = v[:,:,1]
        u[:,:,0] = u[:,:,1]
        b = state.b.view('k')
        b[:, :, 0] = 10.
    if (+1, 0, 0) in ngbs.keys():
        pass
    else:
        w = state.u['k'].view('k')
        w[:, :, -1] = 0.
        v = state.u['j'].view('k')
        u = state.u['i'].view('k')
        v[:,:,-1] = v[:,:,-2]
        u[:,:,-1] = u[:,:,-2]
        b = state.b.view('k')
        b[:, :, -1] = 10.

    if (0, -1, 0) in ngbs.keys():
        pass
    else:
        v = state.u['j'].view('j')
        v[:, :, 0] = 0.
        w = state.u['k'].view('j')
        u = state.u['i'].view('j')
        w[:,:,0] = w[:,:,1]
        u[:,:,0] = u[:,:,1]
        b = state.b.view('j')
        b[:, :, 0] = 10.
    if (0, +1, 0) in ngbs.keys():
        pass
    else:
        v = state.u['j'].view('j')
        v[:, :, -1] = 0.
        w = state.u['k'].view('j')
        u = state.u['i'].view('j')
        w[:,:,-1] = w[:,:,-2]
        u[:,:,-1] = u[:,:,-2]
        b = state.b.view('j')
        b[:, :, -1] = 10.

    if (0, 0, -1) in ngbs.keys():
        pass
    else:
        u = state.u['i'].view('i')
        u[:, :, 0] = 0.
        w = state.u['k'].view('i')
        v = state.u['j'].view('i')
        w[:,:,0] = w[:,:,1]
        v[:,:,0] = v[:,:,1]
        b = state.b.view('i')
        b[:, :, 0] = 10.
    if (0, 0, +1) in ngbs.keys():
        pass
    else:
        u = state.u['i'].view('i')
        u[:, :, -1] = 0.
        w = state.u['k'].view('i')
        v = state.u['j'].view('i')
        w[:,:,-1] = w[:,:,-2]
        v[:,:,-1] = v[:,:,-2]
        b = state.b.view('i')
        b[:, :, -1] = 10.

def apply_bc_on_p(state, ngbs):
    "Pressure continuity on walls"
    if (-1, 0, 0) in ngbs.keys():
        pass
    else:
        p = state.p.view('k')
        p[:,:,0] = p[:,:,1]
    if (+1, 0, 0) in ngbs.keys():
        pass
    else:
        p = state.p.view('k')
        p[:,:,-1] = p[:,:,-2]

    if (0, -1, 0) in ngbs.keys():
        pass
    else:
        p = state.p.view('j')
        p[:, :, 0] = p[:,:,1]
    if (0, +1, 0) in ngbs.keys():
        pass
    else:
        p = state.p.view('j')
        p[:, :, -1] = p[:,:,-2]

    if (0, 0, -1) in ngbs.keys():
        pass
    else:
        p = state.p.view('i')
        p[:, :, 0] = p[:,:,1]

    if (0, 0, +1) in ngbs.keys():
        pass
    else:
        p = state.p.view('i')
        p[:, :, -1] = p[:,:,-2]


def apply_bc_on_vorticity(state, ngbs):
    """ apply free-slip BC on left walls, if any """

    if (-1, 0, 0) in ngbs.keys():
        pass
    else:
        wj = state.vor['j'].view('k')
        wj[:, :, 0] = 0.
        wi = state.vor['i'].view('k')
        wi[:, :, 0] = 0.

    if (0, -1, 0) in ngbs.keys():
        pass
    else:
        wk = state.vor['k'].view('j')
        wk[:, :, 0] = 0.
        wi = state.vor['i'].view('j')
        wi[:, :, 0] = 0.


    if (0, 0, -1) in ngbs.keys():
        pass
    else:
        wk = state.vor['k'].view('i')
        wk[:, :, 0] = 0.
        wj = state.vor['j'].view('i')
        wj[:, :, 0] = 0.
