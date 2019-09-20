"""

The class defining the 3D arrays and the model state

    - Scalars
    - Vectors
    - State

"""
import numpy as np

modelvar = {
    'b': ['scalar', 'buoyancy', 'm.s^-2', 'T'],
    'p': ['scalar', 'pressure', 'm.s^-2', 'T'],
    'ke': ['scalar', 'kinetic energy', 'm^2.s^-2', 'T'],
    'u': ['velocity', 'covariant velocity', 'm^2.s^-1'],
    'U': ['velocity', 'contravariant velocity',  's^-1'],
    'vor': ['vorticity', 'vorticity',  'm^2.s^-1']}

# ----------------------------------------------------------------------


class Scalar(object):
    def __init__(self, param, header):
        """

        header is a dictionnary containing meta information about the scalar
        header.keys = ['name', 'dimensions', 'location', 'nickname']
        this information will be written in the Netcdf file

        """
        self.list_param = ['nx', 'ny', 'nz', 'nh']
        # later 'param' will be an object Param (like in Fluid2d)
        # currently we will suppose that param is a dictionnary
        #
        p = {}
        for k in self.list_param:
            p[k] = param[k]
        self.param = p
#            setattr(self, k, param[k])

        self.header = header
        required_keys = ['name', 'units', 'location', 'nickname']
        keypresent = [k in header.keys() for k in required_keys]
        assert all(keypresent), "header is incomplete"

        self.nature = 'scalar'

        # nx, ny, nz, nh = self.nx, self.ny, self.nz, self.nh
        nx, ny, nz, nh = p['nx'], p['ny'], p['nz'], p['nh']
        self.data = {}
        # extended arrays with halo
        # we might be smarter by removing the halo
        # in the directions where we don't need a halo
        self.data['j'] = np.zeros((nx+2*nh, nz+2*nh, ny+2*nh))
        self.data['k'] = np.zeros((ny+2*nh, nx+2*nh, nz+2*nh))
        self.data['i'] = np.zeros((nz+2*nh, ny+2*nh, nx+2*nh))
        self.activeview = 'i'

    def duplicate(self):
        """ return a new Scalar instance, based on self """
        # p = {'nx': self.nx, 'ny': self.ny, 'nz': self.nz, 'nh': self.nh}
        return Scalar(self.param, self.header)

    def view(self, idx):
        """ return the 3D array

        if lastview is idx
        - copy it from the last activeview buffer to the 'idx' buffer
        - return the pointer to the data

        in any case, the function returns the pointer to the buffer
        """

        if self.activeview == idx:
            field = self.data[idx]

        else:
            current = self.data[self.activeview]
            field = self.data[idx]
            if self.activeview+idx in ['ij', 'jk', 'ki']:
                rightswap(current, field)

            else:
                leftswap(current, field)

            self.activeview = idx

        return field

# ----------------------------------------------------------------------


class Vector(object):
    def __init__(self, param, header, nature='vel'):
        """
        A vector is defined by its three components along the three
        directions i, j and k. We also call these directions x, y and z

        The nature of a vector
        - either at velocity, nature='vel'. It lives at face centers,
        - or vorticity, nature='vor'. It lives along cell edges.

        """
        if nature == 'vel':
            components = {'i': 'U', 'j': 'V', 'k': 'W'}

        elif nature == 'vor':
            components = {'i': 'VW', 'j': 'WU', 'k': 'UV'}

        else:
            raise ValueError('nature of the vector is unknown')

        # define the three components by looping over components
        # self.i, self.j, self.k
        for k in components.keys():
            header.update({'location': components[k], 'comp': k})
            setattr(self, k, Scalar(param, header))

        # use the param from the last component as the param of Vector
        self.param = getattr(self, k).param
        self.header = header
        self.nature = nature

    def duplicate(self):
        return Vector(self.param, self.header, self.nature)

    def view(self, idx):
        """ return the three components """
        return [self.x.view(idx), self.y.view(idx), self.z.view(idx)]

# ----------------------------------------------------------------------


class State(object):
    """
    pack a list of variables (Scalars and Vectors) into a 'state' dictionary

    there is one entry per variable (Scalar) or per components (for Vector)

    the entry key is the variable 'nickname' (Scalar)
    or 'nickname_comp' (comp = 'i', 'j', 'k')

    """

    def __init__(self, listvar):
        # table of content
        toc = {}
        state = {}
        for var in listvar:
            key = var.header['nickname']
            toc[key] = var.nature
            setattr(self, key, var)  # <- store the variable

        self.toc = toc
        self.state = state

    def get(self, var):
        return getattr(self, var)

    def duplicate(self):
        """

        return a new state (with new allocated arrays)

        """
        listvar = []
        for var in self.toc.keys():
            listvar += [getattr(self, var).duplicate()]

        return State(listvar)

    def _print(self):
        for var in self.toc:
            print('%10s : ' % var, getattr(self, var))
        # for k in self.state.keys():
        #     print('%10s : ' % k, self.state[k])

# ----------------------------------------------------------------------


def rightswap(src, dest):
    dest[:, :, :] = np.transpose(src, (2, 0, 1))
    # alternative syntax
    #    dest[:, :, :] = np.moveaxis(src, 2, 0)


def leftswap(src, dest):
    dest[:, :, :] = np.transpose(src, (1, 2, 0))


def get_state(param):
    """
    setup a model state

    useful in the debug phase to rapidly get a full operational state

    """
    listvar = []
    for var, atts in modelvar.items():
        if atts[0] == 'scalar':
            header = {'nickname': var,
                      'name': atts[0],
                      'units': atts[1],
                      'location': atts[2]}
            listvar += [Scalar(param, header)]

        elif atts[0] == 'velocity':
            header = {'nickname': var,
                      'name': atts[0],
                      'units': atts[1]}
            listvar += [Vector(param, header)]

        elif atts[0] == 'vorticity':
            header = {'nickname': var,
                      'name': atts[0],
                      'units': atts[1]}
            listvar += [Vector(param, header, nature='vor')]

    return State(listvar)


def get_work(param):
    """ one may want to quickly get a work array, there you go """
    return Scalar(param, {'name': 'work', 'units': 'any',
                          'location': 'any', 'nickname': 'work'})


# ----------------------------------------------------------------------
if __name__ == '__main__':

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}

    # define the model's state
    s = get_state(param)
    # and a RHS
    ds = s.duplicate()

    print('-'*40)
    print('Model state table of content:')
    print(s.toc)

    print('-'*40)
    print('model state')
    s._print()

    print('-'*40)
    print('model rhs:')
    ds._print()

    print('-'*40)
    assert type(
        s.get('b')).__name__ == 'Scalar', 's.get() does not return a Scalar'
    assert type(
        s.get('u')).__name__ == 'Vector', 's.get() does not return a Vector'
    assert type(
        s.get('u').i).__name__ == 'Scalar', 'vector.i does not return a Scalar'