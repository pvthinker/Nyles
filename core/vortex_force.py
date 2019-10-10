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

    The fortran subroutine imposes 0 velocity at the boudaries.

    TODO : make the function call the fortran subroutine only 3 times.
    """

    positive_triplets = {"ijk","jki","kij"}
    negative_triplets = {"kji","ikj","jik"}


    for triplet in positive_triplets :

        #Using the convention of taking the inner index as the index of the upwinding of vorticity
        covariant_component = rhs.u[triplet[0]].view(triplet[2])
        contravariant_component = state.U[triplet[2]].view(triplet[2])
        vorticity_component = state.vor[triplet[1]].view(triplet[2])

        fortran.vortex_force_calc(contravariant_component, vorticity_component, covariant_component, 1.)

    for triplet in negative_triplets :

        covariant_component = rhs.u[triplet[0]].view(triplet[2])
        contravariant_component = state.U[triplet[2]].view(triplet[2]) #3rd component as inner
        vorticity_component = state.vor[triplet[1]].view(triplet[2])

        fortran.vortex_force_calc(-contravariant_component, vorticity_component, covariant_component, 1.)



if __name__ == '__main__':

    import variables as var

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)
    ds = state.duplicate_prognostic_variables()

    vortex_force(state, ds)
