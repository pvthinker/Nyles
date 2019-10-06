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
from mpi import topology as topo
from collections import namedtuple


# Define the attributes of a variable in the LES model:
#  - type can be 'scalar' or 'velocity' or 'vorticity'
#  - name is the name of the physical quantity
#  - unit is the physical unit of the quantity
#  - prognostic is True for prognostic variables and False for diagnostic variables
ModelVariable = namedtuple('ModelVariable', ['type', 'name', 'unit', 'prognostic'])

# Define the variables in the LES model
modelvar = {
    'b': ModelVariable('scalar', 'buoyancy', 'm.s^-2', prognostic=True),
    'p': ModelVariable('scalar', 'pressure', 'm.s^-2', prognostic=False),
    'ke': ModelVariable('scalar', 'kinetic energy', 'm^2.s^-2', prognostic=False),
    'u': ModelVariable('velocity', 'covariant velocity', 'm^2.s^-1', prognostic=True),
    'U': ModelVariable('velocity', 'contravariant velocity',  's^-1', prognostic=False),
    'vor': ModelVariable('vorticity', 'vorticity',  'm^2.s^-1', prognostic=False),
}

# ----------------------------------------------------------------------


class Scalar(object):
    """A scalar field in 3D space representing a physical quantity.

    The physical quantity is defined by a common name, a short name and
    a physical unit, which will be written in the NetCDF file on save.

    Attributes:
     - param: dictionary describing the properties of the 3D space in
        which the scalar field lives, i.e., its size, the number of
        points in the halo and the neighbours
     - name: a common name describing the physical quantity represented
        by the scalar
     - nickname: a short name or a symbol for the physical quantity
     - unit: the physical unit of the quantity
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
    def __init__(self, param, name, nickname, unit, prognostic: bool):
        """Construct a scalar field in 3D space for a physical quantity.

        The arguments of the constructor have the same role as the
        correspondent attributes of the class.  See in the description
        of the class for more information.
        """
        # Copy the necessary information of param to a local dictionary;
        # this will be replaced by an object of a class Param (like in Fluid2d)
        self.param = {k: param[k] for k in ['nx', 'ny', 'nz', 'nh']}
        self.param['neighbours'] = (
            param['neighbours'] if 'neighbours' in param.keys() else topo.noneighbours()
        )

        self.name = name
        self.nickname = nickname
        self.unit = unit
        self.prognostic = prognostic

        # Calculate size needed in every direction, taking into account
        # the halo in the direction where there is a neighbour
        neighbours = self.param['neighbours']
        self.size = {'i': self.param['nx'], 'j': self.param['ny'], 'k': self.param['nz']}
        for direction in 'ijk':
            # Check for neighbour in plus-direction
            if neighbours[direction+'p'] is not None:
                self.size[direction] += self.param['nh']
            # Check for neighbour in minus-direction
            if neighbours[direction+'m'] is not None:
                self.size[direction] += self.param['nh']
        # Create arrays extended by the halo;
        # it might be smarter to remove the halo
        # in the directions where it is not needed
        self.data = {
            'i': np.zeros((self.size['k'], self.size['j'], self.size['i'])),
            'j': np.zeros((self.size['i'], self.size['k'], self.size['j'])),
            'k': np.zeros((self.size['j'], self.size['i'], self.size['k'])),
        }
        self.activeview = 'i'

    def duplicate(self):
        """Return a new empty Scalar based on self.

        The new scalar field has the same meta-information as self
        (param, name, nickname, unit, prognostic).  Its data arrays are
        initialized with zeros.
        """
        return Scalar(self.param, self.name, self.nickname, self.unit, self.prognostic)

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
                'argument idx of Scalar.flipview must be in ["i","j","k"], not ' + repr(idx)
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
     - param, name, nickname, unit, prognostic: see class Scalar
     - is_velocity: Boolean variable to distinguish velocity vectors
        from vorticity vectors

    Methods:
     - duplicate
     - get_nature
    """
    def __init__(self, param, name, nickname, unit, prognostic: bool, is_velocity):
        """Construct a vector field in 3D space for a physical quantity.

        The arguments of the constructor have the same role as the
        correspondent attributes of the class.  See in the description
        of the class for more information.
        """
        # Create a scalar for each of the three components (self['i'],
        # self['j'], self['k']) in a loop over the three directions
        for direction in 'ijk':
            self[direction] = Scalar(
                param,
                name + ' (' + direction + ')',
                nickname + '_' + direction,
                unit,
                prognostic,
            )

        # Use the attribute param of one of the components as the own
        # attribute param to avoid parsing the argument param again
        self.param = self['i'].param

        self.name = name
        self.nickname = nickname
        self.unit = unit
        self.prognostic = prognostic
        self.is_velocity = is_velocity

    def duplicate(self):
        """Return a new empty Vector based on self.

        The new vector field has the same meta-information as self
        (param, name, nickname, unit, prognostic, is_velocity).  The
        data arrays of its components are initialized with zeros.
        """
        return Vector(self.param, self.name, self.nickname, self.unit,
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
        names of the variables and the values are the nature of the
        variables
     - an object of the class state contains one attribute for each
        element of 'toc' with the key of the dictionary-element as the
        name of the attribute

    Methods:
     - duplicate
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

    def duplicate(self):
        """Return a new state with new allocated arrays."""
        listvar = [getattr(self, var).duplicate() for var in self.toc.keys()]
        return State(listvar)

    def __str__(self):
        return "\n".join(['{:10}: {!r}'.format(var, getattr(self, var)) for var in self.toc])

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
                Scalar(param, var.name, nickname, var.unit, var.prognostic)
            )
        elif var.type == 'velocity':
            listvar.append(
                Vector(param, var.name, nickname, var.unit, var.prognostic, is_velocity=True)
            )
        elif var.type == 'vorticity':
            listvar.append(
                Vector(param, var.name, nickname, var.unit, var.prognostic, is_velocity=False)
            )
    return State(listvar)


def get_work(param):
    """Create quickly an array to array to work with.

    This should not be used apart from debugging purposes.
    """
    return Scalar(param, 'work', 'w', 'any', True)


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
    print(s)

    print('-'*40)
    print('model rhs:')
    print(ds)

    print('-'*40)
    assert type(s.b).__name__ == 'Scalar',\
        'state.b does not return a Scalar'
    assert type(s.u).__name__ == 'Vector',\
        'state.u does not return a Vector'
    assert type(s.u['i']).__name__ == 'Scalar',\
        'vector["i"] does not return a Scalar'

    cov = s.u.view()
    contra = s.U.view()
