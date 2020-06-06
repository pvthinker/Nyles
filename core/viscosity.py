import fortran_dissipation as fortran

def add_viscosity(grid, state, dstate, viscosity):
    for comp in "ijk":
        for direc in "ijk":
            coef = viscosity*grid.ids2[direc]
            cova = state.u[comp].view(direc)
            dudt = dstate.u[comp].view(direc)
            if direc == comp:
                # special case to handle the staggering
                fortran.add_laplacian_uxx(cova, dudt, coef)
            else:
                fortran.add_laplacian(cova, dudt, coef)
