"""
Computes the vortex force term

TODO : Test that is performs what is expected

"""

import fortran_vortex_force as fortran
from timing import timing
import numpy as np

@timing
def vortex_force(state, rhs, order):
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
        And the interpolation of the velocity will need to be in j and i


    The fortran subroutine imposes 0 velocity at the boudaries.


    TODO : make the function call the fortran subroutine only 3 times.


    """

    upw_orders = {1,3,5}

    assert(order in upw_orders)



    # for k, j, i in ["kji", "ikj", "jik"]:
    for k, j, i in ["ikj", "jik", "kji"]:
         #Using the convention of taking the inner index as the index of the upwinding of vorticity

        u_i = rhs.u[i].flipview(j)
        U_k = state.U[k].flipview(j)
        w_j = state.vor[j].flipview(j)
        fortran.vortex_force_calc(U_k, w_j, u_i, 1, order)

        u_k = rhs.u[k].flipview(j)
        U_i = state.U[i].flipview(j)
        uk = flipij(u_k)
        Ui = flipij(U_i)
        wj = flipij(w_j)

        fortran.vortex_force_calc(Ui, wj, uk, -1, order)
        u_k[:] = flipij(uk)
 
 
def flipij(phi):
    return np.transpose(phi, [0, 2, 1])



if __name__ == '__main__':
    import variables as var
    import vorticity
    import topology as topo
    import matplotlib.pyplot as plt
    import numpy as np

    procs = [1, 1, 1]
    topo.topology = "closed"
    myrank = 0

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    nx = 16
    ny = 32
    nz = 64

    param = {'nx': nx, 'ny': ny, 'nz': nz, 'nh': 2, 'neighbours' : neighbours}
    state = var.get_state(param)

    Ui = state.U['i'].view('k')
    Uj = state.U['j'].view('k')
    Uk = state.U['k'].view('k')

    Curr = 4
    Uk[...] = Curr

    wi = state.vor['i'].view('k')
    wj = state.vor['j'].view('k')
    wk = state.vor['k'].view('k')

    x = np.linspace(0,nx,nx)
    y = np.linspace(0,ny,ny)

    X, Y = np.meshgrid(x,y)
    x0 = nx//2
    y0 = ny//2

    x1 = nx//2 + 5
    y1 = ny//2 + 5

    X = np.repeat(X[:,:,np.newaxis], nz, axis = 2)
    Y = np.repeat(Y[:,:,np.newaxis], nz, axis = 2)

    Omega = 2
    wj[...] = Omega * np.exp(-((X-x0)**2/20 + (Y-y0)**2/20))

    ds = state.duplicate_prognostic_variables()

    vortex_force(state, ds, 1)

    ui = ds.u['i'].view()
    uj = ds.u['j'].view()
    uk = ds.u['k'].view()

    plt.figure()
    plt.title('ui')
    cs = plt.pcolor(ui[:,:,10])
    #plt.contour(wj[:,:,10], 5, colors = 'k')
    plt.colorbar(cs)

    plt.figure()
    plt.title('uj')
    cs = plt.pcolor(uj[:,:,10])
    #plt.contour(wj[:,:,10], 5, colors = 'k')
    plt.colorbar(cs)

    plt.figure()
    plt.title('ui')
    cs = plt.pcolor(uk[:,:,10])
    #plt.contour(wj[:,:,10], 5, colors = 'k')
    plt.colorbar(cs)

    plt.figure()
    plt.title('wj')
    plt.pcolor(wj[:,:,10])
    plt.colorbar()


    plt.figure()
    plt.plot(ui[:,8,10], label = 'vf_i')
    plt.plot(wj[:,8,10], label = 'w_j')
    plt.legend()

    plt.show()
