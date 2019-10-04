"""
Just a try while I understand how everything works together

Advection of the buoyancy

vol \partial_t\ b = - \partial_i\ (b vol U^i)

"""

import fortran_upwind as fortran

def buoyancyAdv(state) :
    """
        Advection of the buoyancy, following the model of bernoulli.py

    """

    b = state.b
    U = state.U

    """
        Should be dx, dy, dz to calculate the volume. Not sure from where I can get them ??
    """
    vol = b.nx * b.ny * b.nz

    #preallocate the delta buoyancy, is it the right way ?
    db = b.duplicate()

    for d in 'ijk':
        comp = U[d]

        #upwind(trac, u, dtrac, vol, iflag, l, m, n) what is iflag ??
        fortran.upwind(b.view(d), comp.view(d), db.view(d), vol, b.nx, b.ny, b.nz)

    b += db
