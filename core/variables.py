"""Handling and storage of data in Nyles.

This module contains the essential classes for data handling in Nyles.
The most fundamental type is a scalar field living in 3D space and
representing a physical quantity.  It is modelled by the class Scalar,
which stores the discretization of the field in 3D numpy arrays.
A vector field represents a physical quantity with three components and
is modelled by objects of the class Vector.  Each component of the vector
is a scalar field.  The scalars and vectors of Nyles are managed by an
instance of the class State.

Variables in the LES model can be characterised by the way they appear
in the equations:
 - A prognostic variable is a variable, of which the time derivative
   appears explicitly in the equations, for example buoyancy and
   covariant velocity;
 - A diagnostic variable is a variable, of which the time derivative
   does not appear explicitly in the equations; it cannot be integrated
   in time, so its development must be calculated differently.
"""

import numpy as np
import topology as topo
from collections import namedtuple


# Define the attributes of a variable in the LES model:
#  - type can be 'scalar' or 'velocity' or 'vorticity'
#  - name is the name of the physical quantity
#  - dimension is the physical dimension of the quantity;  this is to
#    associate the correct units with each quantity in the history file
#  - prognostic is True for prognostic variables and False for diagnostic variables
ModelVariable = namedtuple(
    'ModelVariable', ['type', 'name', 'dimension', 'prognostic']
)

# Define the variables in the LES model
modelvar = {
    'b': ModelVariable('scalar', 'buoyancy', 'L.T^-2', prognostic=True),
    'p': ModelVariable('scalar', 'pressure', 'L^2.T^-2', prognostic=False),
    'ke': ModelVariable('scalar', 'kinetic energy', 'L^2.T^-2', prognostic=False),
    'div': ModelVariable('scalar', 'divergence', 'T^-1', prognostic=False),
    'u': ModelVariable('velocity', 'covariant velocity', 'L^2.T^-1', prognostic=True),
    'U': ModelVariable('velocity', 'contravariant velocity',  'T^-1', prognostic=False),
    'vor': ModelVariable('vorticity', 'vorticity',  'L^2.T^-1', prognostic=False),
}

# ----------------------------------------------------------------------


class Scalar(object):
    """A scalar field in 3D space representing a physical quantity.

    The physical quantity is defined by a common name, a short name and
    a dimension. This information is written in the NetCDF file on save.

    The data of the scalar field can be accessed via the "view"-methods,
    which return a 3D array.  For performance reasons, the data is
    internally saved three times in three arrays.

    Attributes:
     - param: dictionary describing the properties of the 3D space in
        which the scalar field lives, i.e., its size, the number of
        points in the halo and the neighbours
     - name: a common name describing the physical quantity represented
        by the scalar
     - nickname: a short name or a symbol for the physical quantity
     - dimension: the physical dimension of the quantity
     - prognostic: a Boolean variable to distinguish prognostic from
        diagnostic variables
     - data: dictionary with keys ['i', 'j', 'k']; each key points to a
        3D array which represents the discretization of the scalar
        field; the array is orientated, such that its key is the inner
        direction; note that the shapes of the arrays can differ;
     - activeview: 'i' or 'j' or 'k', stating which array in 'data'
        contains the up-to-date data
     - size: dictionary with keys ['i', 'j', 'k'] which stores the
       number of grid points in the corresponding direction

    Methods:
     - duplicate
     - view
     - flipview
     - viewlike
     - get_nature

    """

    def __init__(self, param, name, nickname, dimension, prognostic: bool=False):
        """Construct a scalar field in 3D space for a physical quantity.

        The arguments of the constructor have the same role as the
        correspondent attributes of the class.  See in the description
        of the class for more information.
        """
        # Copy the necessary information of param to a local dictionary;
        # this will be replaced by an object of a class Param (like in Fluid2d)
        required_param = ['nx', 'ny', 'nz', 'nh', 'neighbours']
        self.param = {k: param[k] for k in required_param}

        self.name = name
        self.nickname = nickname
        self.dimension = dimension
        self.prognostic = prognostic

        # Calculate size needed in every direction, taking into account
        # the halo in the direction where there is a neighbour
        neighbours = self.param['neighbours']
        p = param
        nx, ny, nz, nh = p['nx'], p['ny'], p['nz'], p['nh']
        shape = [nz, ny, nx]
        self.shape = shape
        size, domainindices = topo.get_variable_shape(shape, neighbours, nh)

        nzl, nyl, nxl = size
        self.size = {'i': nxl, 'j': nyl, 'k': nzl}
        self.domainindices = domainindices

        # define self.mg_idx, the MG array index span
        # MG arrays have the same halo policy
        # than Scalars (halo only if necessary)
        # with the difference that halo width is 1 (nh=1) for MG
        # also MG arrays all use the same (k,j,i) convention
        #
        # if 'field' is a view('i') of a Scalar then
        # field[mg_idx] returns an array the exact same size
        # as a MG array
        k0, k1, j0, j1, i0, i1 = domainindices

        startk, endk = max(0, k0-1), min(nzl, k1+1)
        startj, endj = max(0, j0-1), min(nyl, j1+1)
        starti, endi = max(0, i0-1), min(nxl, i1+1)

        kdx = slice(startk, endk)
        jdx = slice(startj, endj)
        idx = slice(starti, endi)

        self.mg_idx = (kdx, jdx, idx)
        self.mg_idx2 = np.array([startk, endk, startj, endj, starti, endi], dtype=int)

        # Create arrays extended by the halo;
        # it might be smarter to remove the halo
        # in the directions where it is not needed
        self.data = {
            'i': np.zeros((nzl, nyl, nxl)),
            'j': np.zeros((nxl, nzl, nyl)),
            'k': np.zeros((nyl, nxl, nzl)),
        }
        self.activeview = 'i'

    def duplicate(self):
        """Return a new empty Scalar based on self.

        The new scalar field has the same meta-information as self
        (param, name, nickname, dimension, prognostic).  Its data
        arrays are initialized with zeros.
        """
        return Scalar(self.param, self.name, self.nickname, self.dimension, self.prognostic)

    def view(self, idx=None):
        """Return a pointer to the data with idx as inner direction.

        Example:
         - view('i') returns data in the convention (k, j, i)

        Argument:
         - idx: either 'i' or 'j' or 'k' to refer to the 3 directions or
            None to take the current orientation (self.activeview)
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
        """Return pointer to the data array with idx as outer direction.

        Example:
         - view('i') returns data in the convention (i, k, j)

        See 'view' for further information.
        """
        if idx == 'i':
            x = self.view('j')
        elif idx == 'j':
            x = self.view('k')
        elif idx == 'k':
            x = self.view('i')
        else:
            raise ValueError(
                'argument idx of Scalar.flipview must be in ["i","j","k"], not '
                + repr(idx)
            )
        return x

    def viewlike(self, scalar):
        """Return view in the same convention as scalar."""
        return self.view(scalar.activeview)

    @staticmethod
    def get_nature():
        """Return the string 'scalar'.

        This is used to handle scalars and vectors in the same way.
        """
        return 'scalar'

# ----------------------------------------------------------------------


class Vector(dict):
    """A vector field in 3D space representing a physical quantity.

    A vector has three components along the three directions i, j and k,
    also called x, y and z.  Each component is an instance of the class
    Scalar.  A vector can either represent a velocity (living in face
    centres) or the vorticity (living along cell edges).

    The components of the vector are accessed like elements of a
    dictionary, i.e., when v is an instance of the class Vector, its
    component in x-direction is v['i'], in y-direction it is v['j'] and
    in z-direction it is v['k'].

    Attributes:
     - param, name, nickname, dimension, prognostic: see class Scalar
     - is_velocity: Boolean variable to distinguish velocity vectors
        from vorticity vectors

    Methods:
     - duplicate
     - get_nature
    """

    def __init__(self, param, name, nickname, dimension,
                 prognostic: bool=False, is_velocity: bool=True):
        """Construct a vector field in 3D space for a physical quantity.

        The arguments of the constructor have the same role as the
        correspondent attributes of the class.  See in the description
        of the class for more information.
        """
        # Create a scalar for each of the three components (self['i'],
        # self['j'], self['k']) in a loop over the three directions
        dirname = {"i": "x", "j": "y", "k": "z"}
        for direction in 'ijk':
            self[direction] = Scalar(
                param,
                name + " %s-component" % dirname[direction],
                nickname + '_' + direction,
                dimension,
                prognostic,
            )

        # Use the attribute param of one of the components as the own
        # attribute param to avoid parsing the argument param again
        self.param = self['i'].param

        self.name = name
        self.nickname = nickname
        self.dimension = dimension
        self.prognostic = prognostic
        self.is_velocity = is_velocity

    def duplicate(self):
        """Return a new empty Vector based on self.

        The new vector field has the same meta-information as self
        (param, name, nickname, dimension, prognostic, is_velocity).
        The data arrays of its components are initialized with zeros.
        """
        return Vector(self.param, self.name, self.nickname, self.dimension,
                      self.prognostic, self.is_velocity)

    def get_nature(self):
        """Return the string 'velocity' or 'vorticity'.

        This method is useful to handle scalars and vectors in the same
        way.
        """
        if self.is_velocity:
            return 'velocity'
        else:
            return 'vorticity'

# ----------------------------------------------------------------------


class State(object):
    """Container for the variables representing the physical quantities.

    The State contains a table of contents and the variables which
    represent the physical quantities of the model.  These can be Scalar
    fields or Vector fields.

    Attributes:
     - toc: the table of contents is a dictionary; its keys are the
        nicknames of the variables and the values are the nature of
        the variables
     - an object of the class state contains one attribute for each
        element of “toc” with the key of the dictionary-element as
        the name of the attribute

    Methods:
     - duplicate_prognostic_variables
     - get
     - get_prognostic_variables
     - get_prognostic_scalars
    """

    def __init__(self, listvar):
        """Construct a container for the variables of the model.

        Arguments:
         - listvar: list of objects of the classes Scalar or Vector
            to be stored in the State object
        """
        # Create a table of contents (toc) for the State object
        self.toc = {}
        # Store the given variables as attributes and in the toc
        for var in listvar:
            key = var.nickname
            self.toc[key] = var.get_nature()
            setattr(self, key, var)

    def __str__(self):
        return "\n".join(['{:10}: {!r}'.format(var, getattr(self, var)) for var in self.toc])

    def duplicate_prognostic_variables(self):
        """Return a new state with new arrays for prognostic variables."""
        list_of_variables = []
        for variable_name in self.toc:
            variable = getattr(self, variable_name)
            if variable.prognostic:
                list_of_variables.append(variable.duplicate())
        return State(list_of_variables)

    def get(self, variable):
        """Return a pointer to the variable with the given name.

        The given variable name can either be the nickname of a
        (Scalar or Vector) variable or it can be a string like
        "nickname_i" or "nickname_j" or "nickname_k", in which case
        the corresponding component of the Vector "nickname" is
        returned.

        Note:  whenever possible, use the simpler syntax
            state.nickname               # good!
        instead of
            state.get('nickname')        # avoid!
        or, to access a component of a Vector, use
            state.nickname['i']          # good!
        instead of
            state.get('nickname_i')      # avoid!

        """
        # Check if a vector component is requested
        if len(variable) > 2 and variable[-2] == "_" and variable[-1] in 'ijk':
            vector = getattr(self, variable[:-2])
            return vector[variable[-1]]
        else:
            # Otherwise return the requested variable
            return getattr(self, variable)

    def get_prognostic_variables(self):
        """Return a list of the nicknames of the prognostic variables."""
        return [v for v in self.toc if getattr(self, v).prognostic]

    def get_prognostic_scalars(self):
        """Return a list of names of prognostic scalars.

        The list contains:
         - for every prognostic Scalar variable its nickname and
         - for every prognostic Vector variable the strings "nickname_i"
           and "nickname_j" and "nickname_k", referring to its three
           components, with "nickname" replaced by the nickname of the
           Vector variable.

        Each of the strings in this list can be passed as an argument
        to the "get" method to receive a pointer to the Scalar object.

        """
        prognostic_scalars = []
        for nickname in self.get_prognostic_variables():
            if self.toc[nickname] == "scalar":
                prognostic_scalars.append(nickname)
            else:
                prognostic_scalars += [
                    "{}_{}".format(nickname, direction) for direction in "ijk"
                ]
        return prognostic_scalars

# ----------------------------------------------------------------------


def rightswap(src, dest):
    dest[:, :, :] = np.transpose(src, (2, 0, 1))
    # alternative syntax
    #    dest[:, :, :] = np.moveaxis(src, 2, 0)


def leftswap(src, dest):
    dest[:, :, :] = np.transpose(src, (1, 2, 0))


def get_state(param):
    """Setup a model state.

    Useful in the debug phase to rapidly get a full operational state.
    """
    listvar = []
    for nickname, var in modelvar.items():
        if var.type == 'scalar':
            listvar.append(
                Scalar(param, var.name, nickname, var.dimension, var.prognostic)
            )
        elif var.type == 'velocity':
            listvar.append(
                Vector(param, var.name, nickname, var.dimension, var.prognostic,
                       is_velocity=True)
            )
        elif var.type == 'vorticity':
            listvar.append(
                Vector(param, var.name, nickname, var.dimension, var.prognostic,
                       is_velocity=False)
            )
    return State(listvar)


def get_work(param):
    """Create quickly an array to work with.

    This should not be used apart from debugging purposes.
    """
    return Scalar(param, 'work', 'w', 'any', True)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    procs = [4, 2, 1]
    topo.topology = 'closed'
    myrank = 3
    nh = 3

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': nh,
             'neighbours': neighbours}

    # define the model's state
    s = get_state(param)
    # and a RHS
    ds = s.duplicate_prognostic_variables()

    print('-'*40)
    print('Model state table of content:')
    print(s.toc)

    print('-'*40)
    print('model state')
    print(s)

    print('-'*40)
    print('model rhs:')
    print(ds)

    print('-'*40)
    assert type(s.b).__name__ == 'Scalar',\
        'state.b is not a Scalar'
    assert type(s.u).__name__ == 'Vector',\
        'state.u is not a Vector'
    assert type(s.u['i']).__name__ == 'Scalar',\
        'u["i"] is not a Scalar'
    print("no errors in type checking")

    # Check that the list of prognostic scalars is correctly created
    # and that they can be accessed with the get method.
    print('-'*40)
    print("Prognostic variables:")
    for variable in s.get_prognostic_variables():
        print(" - {:3}:".format(variable), s.get(variable))
    print("Prognostic scalars:")
    for scalar in s.get_prognostic_scalars():
        print(" - {:3}:".format(scalar), s.get(scalar))

    b = s.b.view('i')
    print('myrank = %i / loc = ' % myrank, loc, ' / procs = ', procs)
    print('neighbours: ', neighbours)
    print('shape of the full b             : ', np.shape(b))
    print('shape of the b matching MG array: ', np.shape(b[s.b.mg_idx]))
    print('corresponding index span: ', s.b.mg_idx)
