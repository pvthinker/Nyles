"""

The class defining the 3D arrays and the model state

    - Scalars
    - Vectors
    - State

"""
import numpy as np
from mpi import topology as topo

modelvar = {
    'b': ['scalar', 'buoyancy', 'm.s^-2', 'T'],
    'p': ['scalar', 'pressure', 'm.s^-2', 'T'],
    'ke': ['scalar', 'kinetic energy', 'm^2.s^-2', 'T'],
    'u': ['velocity', 'covariant velocity', 'm^2.s^-1'],
    'U': ['velocity', 'contravariant velocity',  's^-1'],
    'vor': ['vorticity', 'vorticity',  'm^2.s^-1'],
}

# ----------------------------------------------------------------------


class Scalar(object):
    def __init__(self, param, header):
        """

        header is a dictionary containing meta information about the
        scalar with the keys ['name', 'units', 'location', 'nickname'].
        This information will be written in the Netcdf file.

        """
        self.list_param = ['nx', 'ny', 'nz', 'nh']
        # later 'param' will be an object Param (like in Fluid2d)
        # currently we will suppose that param is a dictionary
        #
        p = {}
        for k in self.list_param:
            p[k] = param[k]
        key = 'neighbours'
        if key in param.keys():
            p[key] = param[key]
        else:
            p[key] = topo.noneighbours()

        self.param = p
# MR: why is this line commented? Should it be removed?
#            setattr(self, k, param[k])

        # Check header for completeness
        for k in ['name', 'units', 'location', 'nickname']:
            assert k in header.keys(), "header is missing the key " + repr(k)
        self.header = header

        self.nature = 'scalar'

        # nx, ny, nz, nh = self.nx, self.ny, self.nz, self.nh
        nx, ny, nz, nh = p['nx'], p['ny'], p['nz'], p['nh']

        ngbs = p['neighbours']
        size = [nx, ny, nz]
        for k, direc in enumerate('ijk'):
            for pm in 'mp':
                if ngbs[direc+pm] is None:
                    pass
                else:
                    size[k] += nh
        nxl, nyl, nzl = size
        self.size = size
        self.data = {}
        # extended arrays with halo
        # we might be smarter by removing the halo
        # in the directions where we don't need a halo
        self.data['j'] = np.zeros((nxl, nzl, nyl))
        self.data['k'] = np.zeros((nyl, nxl, nzl))
        self.data['i'] = np.zeros((nzl, nyl, nxl))
        self.activeview = 'i'

    def duplicate(self):
        """ return a new Scalar instance, based on self """
        # p = {'nx': self.nx, 'ny': self.ny, 'nz': self.nz, 'nh': self.nh}
        return Scalar(self.param, self.header)

    def view(self, idx=None):
        """ return the 3D array with idx as inner direction

        idx = 'i', 'j', 'k'

        view('i') return x in the convention (k, j, i)

        if lastview is idx
        - copy it from the last activeview buffer to the 'idx' buffer
        - return the pointer to the data

        in any case, the function returns the pointer to the buffer

        """
        if idx == self.activeview or idx is None:
            field = self.data[self.activeview]

        else:
            current = self.data[self.activeview]
            field = self.data[idx]
            if self.activeview+idx in ['ij', 'jk', 'ki']:
                rightswap(current, field)

            else:
                leftswap(current, field)

            self.activeview = idx

        return field

    def flipview(self, idx):
        """

        return the 3D array but with idx indicating the outer direction

        flipview('i') return x in the convention (i, k, j)

        """
        if idx == 'i':
            x = self.view('j')
        if idx == 'j':
            x = self.view('k')
        if idx == 'k':
            x = self.view('i')
        return x

    def viewlike(self, scalar):
        """
        return the view that has the same convention than scalar
        """
        activeview = scalar.activeview
        return self.view(activeview)

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
        # MR: This line seems to depend on the order in which the keys
        #     are given to the dictionary.  If that is the case, it
        #     should NOT be done like this, because until Python 3.6,
        #     the insertion order of a dictionary is not preserved.
        #     Suggestion: if the last element is always "k", just write
        #     self.k.param -- otherwise, use a different datastructure
        #     like an OrderedDict or create an extra variable for the
        #     last component.
        #     If this line does not depend on the insertion order, then
        #     just do self.k.param, since it is much clearer.
        self.param = getattr(self, k).param
        self.header = header
        self.nature = nature

    def duplicate(self):
        return Vector(self.param, self.header, self.nature)

    def view(self, idx=None):
        """ return the three components """
        if idx is None:
            idx = self.i.activeview
        return [self.i.view(idx), self.j.view(idx), self.k.view(idx)]

# ----------------------------------------------------------------------


class State(object):
    """
    Pack a list of variables (Scalars and Vectors) into a 'state' dictionary

    The dictionary 'state' has
     - for scalars: one entry per scalar with the key 'nickname'
     - for vectors: one entry per component with the key 'nickname_comp'
       with 'comp' = 'i', 'j' or 'k'

    MR: the information given in this docstring does not represent the
        way it is currently implemented.
    """

    def __init__(self, listvar):
        # Create a table of contents (toc) for the State object
        self.toc = {}
        for var in listvar:
            key = var.header['nickname']
            self.toc[key] = var.nature
            setattr(self, key, var)  # <- store the variable
        self.state = {}

    def get(self, var):
        # Check if 'var' is of the scheme 'nickname_i' (or with '_j' or '_k').
        if len(var) > 2 and var[-2] == "_" and var[-1] in 'ijk':
            # If so, take the vector 'nickname' and return its component 'i'.
            vector = getattr(self, var[:-2])
            return getattr(vector, var[-1])
        else:
            # Otherwise, return the requested variable.
            return getattr(self, var)

    def duplicate(self):
        """

        return a new state (with new allocated arrays)

        """
        listvar = [getattr(self, var).duplicate() for var in self.toc.keys()]

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

    procs = [4, 2, 1]
    topology = 'closed'
    myrank = 3
    nh = 3

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs, topology)

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': nh,
             'neighbours': neighbours}

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
    assert type(
        s.get('u_i')).__name__ == 'Scalar', 'vector_i does not return a Scalar'

    cov = s.get('u').view()
    contra = s.get('U').view()
