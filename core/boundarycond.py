import numpy as np

def apply_bc_on_velocity(state, ngbs):
    """ apply no-flow BC on left walls, if any """
    if (-1, 0, 0) in ngbs.keys():
        pass
    else:
        w = state.u['k'].view('k')
        w[:, :, 0] = 0.
        b = state.b.view('k')
        b[:, :, 0] = 99.

    if (0, -1, 0) in ngbs.keys():
        pass
    else: 
        v = state.u['j'].view('j')
        v[:, :, 0] = 0.
        b = state.b.view('j')
        b[:, :, 0] = 99.
        
    if (0, 0, -1) in ngbs.keys():
        pass
    else:
        u = state.u['i'].view('i')
        u[:, :, 0] = 0.
        b = state.b.view('i')
        b[:, :, 0] = 99.

        
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
