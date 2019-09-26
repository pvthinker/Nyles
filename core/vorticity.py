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

    cova = state.get('u')
    vort = state.get('vor')

    perm = {'i': ('j', 'k'),
            'j': ('k', 'i'),
            'k': ('i', 'j')}

    for diri in 'ijk':
        dirj, dirk = perm[diri]

        uj = getattr(cova, dirj).view(dirj)
        uk = getattr(cova, dirk).view(dirj)
        wi = getattr(vort, diri).view(dirj)

        fortran.vorticity(uj, uk, wi)


@timing
def vorticity_all_comp(state):
    """
    purpose : reference point to compare with the previous

    verdict: it's faster!
    """

    cova = state.get('u')
    vort = state.get('vor')

    ui = cova.i.view('i')
    uj = cova.j.view('i')
    uk = cova.k.view('i')

    wi = vort.i.view('i')
    wj = vort.j.view('i')
    wk = vort.k.view('i')

    fortran.vorticity_all_comp(ui, uj, uk, wi, wj, wk)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)

    vorticity(state)

    vorticity_all_comp(state)
