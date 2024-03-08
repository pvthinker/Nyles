import fortran_dissipation as fortran

def add_viscosity(grid, state, dstate, viscosity):
    for direc in "ijk":
        coef = viscosity*grid.ids2[direc]
        for comp in "ijk":
            cova = state.u[comp].view(direc)
            dudt = dstate.u[comp].view(direc)
            fortran.add_laplacian(cova, dudt, coef)
