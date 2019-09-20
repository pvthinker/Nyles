import fortran_kinenergy as fortran
import variables as var

def kinenergy(state):
    """ 
    compute the kinetic energy function from the model state

    """
    # loop over direction
    for k, direc in enumerate('ijk'):
        u = state.state['u_'+direc].view(direc)
        
        # in sigma coordinates, use U, the contravariant velocity
        # U = state.state['U_'+direc].view(direc)
        
        ke = state.state['ke'].view(direc)
        # in z coordinates, the kinetic energy is ke = u**2 / ds**2
        # where ds2 is the diagonal term of the inverse metric tensor
        ds2 = 1. # TODO: use dx**-2, dy**-2, dz**-2
        
        if k == 0:
            fortran.kin(u, u, ke, ds2, 1) # overwrite ke
            
        else:
            fortran.kin(u, u, ke, ds2, 0) # add to ke
            
            
#----------------------------------------------------------------------
if __name__ == '__main__':

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)

    kinenergy(state)
    
