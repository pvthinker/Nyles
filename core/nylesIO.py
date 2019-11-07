"""Provide Input/Output handling for Nyles through the class NylesIO.

The class NylesIO serves for two closely related purposes:
 1. Output: Store the output of a model run in a netCDF file.
 2. Input:  Load the output of a previously saved model run and continue
            the experiment from the last saved point.

For every experiment, a new directory is created, which contains all the
output files for this model run, referred to as the "output directory".
All the output directories are located in one directory, referred to as
the "data directory".  Currently, the only output file is a netCDF
"history file", which contains screenshots of 3D scalar variables of the
model and an overview of all the experiment parameters.  The experiment
parameters may be shortened if they are too long and not of any of the
following types: bool, int, float, str.

The name of the output directory is given by the name of the experiment,
possibly extended by a counter.  The same name is used for the history
file, extended by "_hist.nc".

The NylesIO class supports three modes:
 - overwrite (default): the experiment name is used as it is, without
    checking if a output directory with the same name exists already
 - count: add a counter to the end of the experiment name to guarantee
    uniqueness of the output directory
 - continue: load a previously saved model run and continue the
    experiment from the last saved point; the experiment name is
    extended by a continue-counter, or this counter is increased if it
    already exists; (not yet implemented)
    Caution: continuing the experiment with the continue-counter 1 will
    take the counter value 2, even if such an output directory exists
    already; its content will be overwritten!  This is to ensure that a
    higher counter is always the continuation of a lower counter.

TODO:
 - add support for multi-core runs
 - complete support for continuing experiments
 - implement diagnostic and flux file if necessary
 - save the array without its halo
"""

import os
import netCDF4 as nc
from enum import Enum   # since Python 3.4
from collections import namedtuple
# Local import
from variables import State


# Define the attributes of a netCDF variable
NetCDFVariable = namedtuple('NetCDFVariable', ['nickname', 'name', 'unit'])


class NylesIO(object):
    """Input/Output handler for Nyles.

    The interaction with an instance of this class should usually only
    occur via the public attributes and methods listed below, and not by
    accessing the private attributes or methods.

    Attributes for public access:
     - output_directory: path to the output directory
     - hist_path: path to the history file; this string contains the
        path to the output directory

    Attributes mostly for private access:
     - hist_variables: list of NetCDFVariable instances, describing the
        scalar variables of the model to be saved in the history file;
        this list is created in the constructor, but its content is
        neither checked nor finalised until "init" is called
     - experiment_parameters: dictionary with all parameters of the
        experiment to be saved in the history file, which are of type
        int or float or string
     - dt_hist: time between two savings in the history file
     - t_next_hist: point in time when to create the next entry in the
        history file
     - n_hist: counter for the number of entries in the history file;
        this value can be readout to display a status message

    Constants mostly for private access:
     - MAX_LENGTH_ATTRIBUTE: maximal number of characters for attributes
        saved in the history file; this applies only to experiment
        parameters which are not of any of the following types:
        bool, int, float, str

    Methods for public access:
     - init
     - do
     - save_array_3D

    Methods mostly for private access:
     - create_history_file
     - write_history_file
    """

    MAX_LENGTH_ATTRIBUTE = 100

    def __init__(self, param):
        """Create an object to handle Input and Output.

        The argument "param" determines the behaviour of this object and
        provides the experiment parameters which are saved in the
        history file.

        The following keys must be present in "param":
         ~ variables_in_history: list of nicknames of the scalar
            variables to be saved in the history file or "all" to save
            all model variables except "work"
         - timestep_history: target timestep for savings in history file
         - datadir: path to the data directory
         - expname: name of the experiment

        The following key can be present in "param":
         - mode: can be "overwrite" or "count" or "continue"; if the key
            is not present, "overwrite" is used; "continue" is not yet
            implemented

        Every key-value-pair in "param" is saved in the history file.
        Values of type int, float or string are saved as they are.
        Boolean variables are converted to the string "True" or "False".
        Every other variable is converted to a string and shortened to
        MAX_LENGTH_ATTRIBUTE characters if necessary.
        """
        self.hist_variables = param["variables_in_history"]  # this will be finalised on "init"
        self.dt_hist = param["timestep_history"]
        self.t_next_hist = 0.0
        self.n_hist = 0

        # Create a copy of the experiment parameters to save them it in the history file
        # TODO: simplify the next statement, if the implementation is settled
        # on param being a dictionary or param being an object of a Param class.
        self.experiment_parameters = param.copy() if isinstance(param, dict) else vars(param)
        # Convert experiment parameters to something that can be saved in a netCDF file
        for key, value in self.experiment_parameters.items():
            if isinstance(value, bool):
                # Convert True and False to strings for saving.
                # Caution: you cannot use bool("False") to convert them back to Boolean type!
                self.experiment_parameters[key] = str(value)
            elif isinstance(value, (int, float, str)):
                # These types can be saved in the netCDF file directly.
                continue
            else:
                # Convert parameters of any other type to a string
                # denoting its type and a representation of its value
                string_representation = (
                    "{}: {}".format(type(value), repr(value))[:self.MAX_LENGTH_ATTRIBUTE+1]
                )
                # Cut strings which are too long
                if len(string_representation) > self.MAX_LENGTH_ATTRIBUTE:
                    string_representation = (
                        string_representation[:self.MAX_LENGTH_ATTRIBUTE-3] + "..."
                    )
                self.experiment_parameters[key] = string_representation

        # Create paths for the in- and output
        datadir = os.path.expanduser(param["datadir"])
        expname = param["expname"]
        if "/" in expname:
            raise ValueError('expname may not contain a slash: "/"')
        out_dir = os.path.join(datadir, expname)
        if "mode" not in param or param["mode"] == "overwrite":
            # Nothing to do, take the path as it is
            pass
        elif param["mode"] == "count":
            counter = 0
            full_expname = "{}_{:04d}".format(expname, counter)
            out_dir = os.path.join(datadir, full_expname)
            while os.path.exists(out_dir):
                counter += 1
                full_expname = "{}_{:04d}".format(expname, counter)
                out_dir = os.path.join(datadir, full_expname)
            expname = full_expname
        elif param["mode"] == "continue":
            if not os.path.exists(out_dir):
                raise ValueError(
                    "cannot load experiment to continue: no directory exists "
                    "with the name {!r}; please specify a correct path or change "
                    "the mode of the experiment.".format(out_dir)
                )
            # TODO
            raise NotImplementedError
        else:
            raise ValueError("unknown mode: {}".format(param["mode"]))
        self.hist_path = os.path.join(out_dir, expname + "_hist.nc")
        self.output_directory = out_dir

    def init(self, state: State, t=0.0, n=0):
        """Initialise the I/O object and create the history file.

        The list of variables to be saved in the history file is
        finalised by either checking if every given variable exists in
        the state or by replacing "all" by a list of all the variables
        in the state (without the variable "work" if it exists).  Every
        vector is split up into its components, so that the list
        contains only scalar variables.

        A new history file is created and the initial state of the model
        is saved.

        Arguments:
         - state: current (initial) state of the model run
         - t: initial time of the model run, usually 0.0, but can be
            higher if a previous experiment is continued
         - n: integration step of the model, usually 0
        """
        # Finalise the list of variables to save in the history file
        if self.hist_variables == "all":
            self.hist_variables = [v for v in state.toc.keys() if v != "work"]
        # Check each variable for correctness and split vectors into their components
        full_hist_variables = []
        for variable in self.hist_variables:
            if variable not in state.toc:
                raise ValueError(
                    "unknown variable {!r} in parameter 'variables_in_history'."
                    .format(variable)
                )
            elif state.toc[variable] == "scalar":
                v = state.get(variable)
                full_hist_variables.append(
                    NetCDFVariable(variable, v.name, v.unit)
                )
            else:
                for component in "ijk":
                    nickname = "{}_{}".format(variable, component)
                    v = state.get(nickname)
                    full_hist_variables.append(
                        NetCDFVariable(nickname, v.name, v.unit)
                    )
        self.hist_variables = full_hist_variables[:]
        # Print a warning if there are no variables to save, but do
        # not stop, because this could be the desired behaviour
        if not self.hist_variables:
            print(
                "Warning: no 3D data will be saved in the history file! "
                "Modify the parameter 'variables_in_history' to change this."
            )
        # Create the output directory if necessary
        if not os.path.isdir(self.output_directory):
            os.makedirs(self.output_directory)
        # Take any model variable to determine the size of any data array
        dimensions = state.get(list(state.toc.keys())[0]).size
        self.create_history_file(dimensions)
        self.write_history_file(state, t, n)
        self.t_next_hist = t + self.dt_hist

    def do(self, state, t, n):
        """Write the current state to the history file if it is time.

        Arguments:
         - state: current state of the model run
         - t: current time of the model run
         - n: current integration step of the model
        """
        # Write data to the history file if it is time to do so
        if t >= self.t_next_hist:
            self.write_history_file(state, t, n)
            self.t_next_hist += self.dt_hist

    def save_array_3D(self, data, name, description=""):
        """Save a 3D array of data in the history file.

        Arguments:
         - data: 3D array of floats to be saved in the history file;
            it must be given in the convention [z, y, x], which is
            returned by view("i") for variables of type Scalar
         - name: a short name used to reference the variable
         - description (optional): additional text saved in the history
            file together with the variable
        """
        with nc.Dataset(self.hist_path, "a") as ncfile:
            v = ncfile.createVariable(name, float, ("z", "y", "x"))
            if description:
                v.long_name = description
            ncfile[name] = data

    def create_history_file(self, dimensions):
        """Create a new history file.

        The created history file contains six variables (x, y, z, t, n)
        with one dimension each, and the variables marked for saving in
        the history file with four dimensions (t, z, y, x) each.  Every
        variable must have a name and can have a description and can
        have a unit.

        The created history file contains all the experiment parameters
        as global attributes.

        Attributes:
         - dimensions: dictionary with keys ['i', 'j', 'k'] and as
            values the number of grid points in the corresponding
            direction
        """
        # Create a new history file for writing
        with nc.Dataset(self.hist_path, "w") as ncfile:
            # Store the experiment parameters
            ncfile.setncatts(self.experiment_parameters)

            # Create the dimensions
            ncfile.createDimension("t")  # unlimited size
            ncfile.createDimension("x", dimensions["i"])
            ncfile.createDimension("y", dimensions["j"])
            ncfile.createDimension("z", dimensions["k"])

            # Create the variables with one dimension
            v = ncfile.createVariable("n", int, ("t",))
            v.long_name = "integration step in the model run"

            # TODO: make the units of x,y,z,t more general
            v = ncfile.createVariable("t", float, ("t",))
            v.long_name = "time in the model run"
            v.units = "s"

            # TODO: explain where the coordinates are taken (cell centre or other)
            v = ncfile.createVariable("x", float, ("x",))
            v.long_name = "coordinate in the x-direction"
            v.units = "m"

            v = ncfile.createVariable("y", float, ("y",))
            v.long_name = "coordinate in the y-direction"
            v.units = "m"

            v = ncfile.createVariable("z", float, ("z",))
            v.long_name = "coordinate in the z-direction"
            v.units = "m"

            # TODO: add mask if a mask is implemented

            # Create variables for the model data.
            # For more information on the attributes of netCDF variables, see
            # https://www.unidata.ucar.edu/software/netcdf/docs/attribute_conventions.html
            for variable in self.hist_variables:
                # Spatial dimensions are in reversed order, because
                # the arrays are stored like this in the Scalar class
                v = ncfile.createVariable(variable.nickname, float, ("t", "z", "y", "x"))
                v.long_name = variable.name
                v.units = variable.unit

            # Save the coordinates
            # TODO: implement this; I currently don't know how to access the
            # grid in Nyles, so here is what it looked like in fluid2d:
            # ncfile['x'][:] = self.grid.xr[0, nh:-nh]
            # ncfile['y'][:] = self.grid.yr[nh:-nh, 0]

    def write_history_file(self, state, t, n):
        """Append the given state to the history file.

        Arguments:
         - state: state of the model run to be saved in the history file
         - t: time of the model run at which this state is taken
         - n: integration step of the model at which this state is taken
        """
        # Open the history file for appending data to it
        with nc.Dataset(self.hist_path, "a") as ncfile:
            ncfile["t"][self.n_hist] = t
            ncfile["n"][self.n_hist] = n
            for v in self.hist_variables:
                ncfile[v.nickname][self.n_hist] = state.get(v.nickname).view("i")
        self.n_hist += 1


if __name__ == "__main__":
    from variables import get_state
    import topology as topo
    import numpy as np

    procs = [4, 2, 1]
    topo.topology = 'closed'
    myrank = 3

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    param = {
        ### Parameters necessary for the IO class
        "datadir": "~/data/Nyles",
        "expname": "test_exp",
        "timestep_history": 1.0,
        ## Choose between a list or "all" ("all" does not include "work")
        "variables_in_history": "all",
        # "variables_in_history": ["u", "b"],
        ## Select one or zero of the following options
        "mode": "overwrite",
        # "mode": "count",
        # "mode": "continue",
        ### Parameters necessary for State
        "nx": 5,
        "ny": 4,
        "nz": 2,
        "nh": 2,
        "neighbours": neighbours,
        ### Parameters to test the storage of parameters in the history file
        "a boolean variable": True,
        "a 2nd boolean variable": False,
        "a long list": list(range(50)),
    }
    state = get_state(param)
    print("* Full state:")
    print(state)

    io = NylesIO(param)
    print("* Output directory:", io.output_directory)
    print("* History file:", io.hist_path)
    print("* Experiment parameters to save:")
    print(io.experiment_parameters)
    io.init(state)
    print("* Variables in the history file:")
    print(io.hist_variables)
    b = state.get("b").view("i")
    for i, t in enumerate([0.2, 0.5, 0.7, 1.0, 1.5, 2.1, 3.0]):
        n = i+1
        print("* t = {:9.2f} / n = {:6} / n_hist = {:6}".format(t, n, io.n_hist))
        b[:] = t * np.ones_like(b)
        io.do(state, t, n)
