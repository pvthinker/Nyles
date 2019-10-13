"""
Computes the vortex force term

TODO : Test that is performs what is expected

"""

import fortran_vortex_force as fortran
from timing import timing

@timing
def vortex_force(state, rhs):
    """
    Calculates the vortex force term

    This function uses the known non-zero triplets of the Levi-Civita tensor to make the
    cross product between the vorticity and contravariant velocity. That way, it makes
    6 calls the fortran subroutine (3 positive and 3 negative).

    The inner index of all the arrays is taken to be the index of the direction of upwinding,
    this should make the computation more efficient in fortran.

    for example :
        d/dt (u_i) = - ( eps_ijk * w_j * U^k + eps_ikj * w_k * U^j )

        Lets consider only the first term : eps_ijk * w_j * U^k

        U^k  is in the k direction, so that the upwinding of vorticity will be done also
        in the k direction.
        This way, we need the k direction to be the inner index of the array, they will
        then be in the form :
        U(j,i,k)
        vort(j,i,k)
        u(j,i,k)
        And the interpolation of the velocity will need to be in k and i

        In the case of the second term, the upwinding will be in U^j, so the inner index
        will be the j direction. The arrays will then be in the form :
        U(i,k,j)
        vort(i,k,j)
        u(i,k,j)
        And the interpolation of the velocity will need to be in k and i


    The fortran subroutine imposes 0 velocity at the boudaries.


    TODO : make the function call the fortran subroutine only 3 times.
    """

    positive_triplets = ["ijk","jki","kij"]
    negative_triplets = ["kji","ikj","jik"]


    for i,j,k in positive_triplets :

        #Using the convention of taking the inner index as the index of the upwinding of vorticity
        u_comp = rhs.u[i].view(k)
        U_comp = state.U[k].view(k)
        w_comp = state.vor[j].view(k)

        fortran.vortex_force_calc(U_comp, w_comp, u_comp, 1, 1.)

    for i,j,k in negative_triplets :

        u_comp = rhs.u[i].view(k)
        U_comp = state.U[k].view(k)
        w_comp = state.vor[j].view(k)

        fortran.vortex_force_calc(U_comp, w_comp, u_comp, -1, 1.)



if __name__ == '__main__':

    import variables as var
    import vorticity
    import matplotlib.pyplot as plt
    import numpy as np

    #Verify that the cross product works

    nx = 40
    ny = 50
    nz = 60
    nh = 2

    param = {'nx': nx, 'ny': ny, 'nz': nz, 'nh': nh}
    state = var.get_state(param)

    Ui = state.U['i'].view()
    Uj = state.U['j'].view()
    Uk = state.U['k'].view()

    Uk[:,:,:] = -1

    wi = state.vor['i'].view()
    wj = state.vor['j'].view()
    wk = state.vor['k'].view()

    wj[:,:,:] = 2

    ds = state.duplicate_prognostic_variables()

    vortex_force(state, ds)

    ui = ds.u['i'].view()
    uj = ds.u['j'].view()
    uk = ds.u['k'].view()


    plt.figure()
    plt.title('ui')
    plt.pcolor(ui[10,:,:])
    plt.colorbar()

    plt.figure()
    plt.title('uj')
    plt.pcolor(uj[10,:,:])
    plt.colorbar()

    plt.figure()
    plt.title('uk')
    plt.pcolor(uk[10,:,:])
    plt.colorbar()

    plt.show()
