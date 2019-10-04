import fortran_vorticity as fortran
import variables as var
from timing import timing


@timing
def vorticity(state):
    """
    compute the vorticity

    omega_i = delta_j[u_k] - delta_k[u_j]

    direction i should be the first one,
    imposing to view arrays in the j direction

    TODO:
       add the Coriolis term
       carefully handle the boundary condition: no-slip or free-slip
    """

    perm = {'i': ('j', 'k'),
            'j': ('k', 'i'),
            'k': ('i', 'j')}

    for diri in 'ijk':
        dirj, dirk = perm[diri]

        uj = state.u[dirj].view(dirj)
        uk = state.u[dirk].view(dirj)
        wi = state.vor[diri].view(dirj)

        fortran.vorticity(uj, uk, wi)


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

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)

    vorticity(state)

    vorticity_all_comp(state)
