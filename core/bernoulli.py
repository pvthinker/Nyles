"""

add b*grad(z) + grad(ke) to rhs

TODO:
  to increase computational intensity, compute delta[b*grad(z) + grad(ke)]
  and add it to the rhs of the multigrid

"""

import fortran_bernoulli as fortran
from timing import timing


@timing
def bernoulli(state, rhs, grid):
    for direction in 'ijk':
        u_component = rhs.u[direction]
        if direction in 'ij':
            fortran.gradke(state.ke.view(direction), u_component.view(direction))
        else:
            fortran.gradkeandb(state.ke.view(direction), state.b.view(direction),
                               u_component.view(direction), grid.dz)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    import variables as var
    from grid import Grid


    param = {
        'nx': 40, 'ny': 50, 'nz': 60,
        'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,
        'neighbours': {}, 'nh': 3,
    }
    grid = Grid(param)
    state = var.get_state(param)
    ds = state.duplicate_prognostic_variables()

    bernoulli(state, ds, grid)
