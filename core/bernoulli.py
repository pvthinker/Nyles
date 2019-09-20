"""

add b*grad(z) + grad(ke) to rhs

"""

import fortran_bernoulli as fortran
import variables as var
from timing import timing


@timing
def bernoulli(state, rhs):

    du = rhs.get('u')
    b = state.get('b')
    ke = state.get('ke')

    for d in 'ijk':
        comp = getattr(du, d)
        if d in 'ij':
            fortran.gradke(ke.view(d), comp.view(d))

        else:
            fortran.gradkeandb(ke.view(d), b.view(d), comp.view(d))


# ----------------------------------------------------------------------
if __name__ == '__main__':

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)
    ds = state.duplicate()

    bernoulli(state, ds)
