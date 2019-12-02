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

For more information on the attributes of netCDF variables, visit
https://www.unidata.ucar.edu/software/netcdf/docs/attribute_conventions.html

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
from grid import Grid

# Define the attributes of a netCDF variable
NetCDFVariable = namedtuple('NetCDFVariable', ['nickname', 'name', 'dimension'])


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
     - last_saved_frame: integration step of the last entry that was
        saved in the history file
     - disk_limit: if the available disk space drops below this value,
        display a warning

    Constants mostly for private access:
     - MAX_LENGTH_ATTRIBUTE: maximal number of characters for attributes
        saved in the history file; this applies only to experiment
        parameters which are not of any of the following types:
        bool, int, float, str

    Methods for public access:
     - init
     - do
     - finalize
     - save_array_3D

    Methods mostly for private access:
     - create_history_file
     - write_history_file
     - get_disk_space_in_GB
    """

    MAX_LENGTH_ATTRIBUTE = 100

    def __init__(self, param):
        """Create an object to handle Input and Output.

        The argument "param" determines the behaviour of this object and
        provides the experiment parameters which are saved in the
        history file.

        The following keys must be present in "param":
         ~ variables_in_history: list of nicknames of the variables to
            be saved in the history file; instead of a list, this
            parameter can have the value "all" to save all model
            variables; for vectors, it is discouraged to include
            individual vector components in this list, because they will
            not be saved at the correct coordinates of the staggered
            grid; instead, always use the whole vector "u", "U" or "vor"
         - timestep_history: target timestep for saving the model state
            in the history file; use 0.0 to save every frame;
         - datadir: path to the data directory
         - expname: name of the experiment
         - mode: one out of ["overwrite", "count", "continue"];
            "continue" is not yet implemented
         - disk_space_warning: limit in GB at which to display a warning
            about low disk space
         - unit_length, unit_duration: the physical units for length and
            duration used in the model; to be saved in the history file

        Every key-value-pair in "param" is saved in the history file.
        Values of type int, float or string are saved as they are.
        Boolean variables are converted to the string "True" or "False".
        Every other variable is converted to a string and shortened to
        MAX_LENGTH_ATTRIBUTE characters if necessary.
        """
        self.disk_limit = param["disk_space_warning"]
        self.hist_variables = param["variables_in_history"]  # this will be finalised on "init"
        self.dt_hist = param["timestep_history"]
        self.t_next_hist = 0.0
        self.n_hist = 0
        self.last_saved_frame = None
        # Define a function to replace dimensions by units for convenience
        self.unit = lambda dimension: (
            dimension.replace("T", param["unit_duration"]).replace("L", param["unit_length"])
        )

        # Create a copy of the experiment parameters to save them it in the history file
        self.experiment_parameters = param.copy()
        # Convert experiment parameters to something that can be saved in a netCDF file
        for key, value in self.experiment_parameters.items():
            if isinstance(value, bool):
                # Convert True and False to strings for saving.
                # Caution: it is not possible to use bool() to convert
                # the strings back to Boolean type, because
                # bool("False") returns True
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
        out_dir = os.path.join(datadir, expname)
        if param["mode"] == "overwrite":
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

    def init(self, state: State, grid: Grid, t=0.0, n=0):
        """Initialise the I/O object and create the history file.

        The list of variables to be saved in the history file is
        finalised by either checking if every given variable exists in
        the state or by replacing "all" by a list of all the variables
        in the state.  Every vector is split up into its components, so
        that the list contains only scalar variables.

        A new history file is created and the initial state of the model
        is saved.

        Arguments:
         - state: current (initial) state of the model run
         - t: initial time of the model run, usually 0.0, but can be
            higher if a previous experiment is continued
         - n: integration step of the model, usually 0
        """
        # Check the list of variables to save in the history file
        if not self.hist_variables:
            # Do not throw an error but only print a warning, because it
            # could be the desired to save no 3D fields
            print(
                "Warning: no 3D data will be saved in the history file! "
                "Modify the parameter 'variables_in_history' to change this."
            )
        elif self.hist_variables == "all":
            self.hist_variables = state.toc.keys()
        elif any(variable not in state.toc for variable in self.hist_variables):
            # Throw an error because this is likely due to a typo of the user
            raise ValueError(
                "a variable chosen for the history file does not exist in this "
                "model.  Available variables: " + str(state.toc)
            )
        # Make list of strings into list of Scalar and Vector objects
        self.hist_variables = [state.get(v) for v in self.hist_variables]
        # Create (if necessary) the output directory and the history file
        if not os.path.isdir(self.output_directory):
            os.makedirs(self.output_directory)
        self.create_history_file(grid)
        # Split vectors into their components to make saving data easier
        full_hist_variables = []
        for v in self.hist_variables:
            if v.get_nature() == "scalar":
                full_hist_variables.append(
                    NetCDFVariable(v.nickname, v.name, v.dimension)
                )
            else:
                for i in "ijk":
                    full_hist_variables.append(
                        NetCDFVariable(v[i].nickname, v[i].name, v[i].dimension)
                    )
        self.hist_variables = full_hist_variables[:]
        self.write_history_file(state, t, n)
        self.t_next_hist = t + self.dt_hist

    def do(self, state, t, n):
        """Write the current state to the history file if it is time.

        When new data was written to the disk, check if the available
        disk space is still above the limit.  Otherwise display a
        warning and pause until the user has decided to continue.

        Return True if the user decides to stop the simulation after a
        low-disk-space warning.  Return False if everything is fine.

        Arguments:
         - state: current state of the model run
         - t: current time of the model run
         - n: current integration step of the model
        """
        # Write data to the history file if it is time to do so
        if t >= self.t_next_hist:
            self.write_history_file(state, t, n)
            self.t_next_hist += self.dt_hist
            # Check if the remaining disk space is sufficient
            if self.disk_limit > 0:
                try:
                    free_space = self.get_disk_space_in_GB()
                except Exception as e:
                    print("Cannot determine available disk space.  Error message:", e)
                    print("Disabling further disk space checks.")
                    self.disk_limit = 0
                else:
                    if free_space < self.disk_limit:
                        print("!"*50)
                        print("Warning, low disk space:")
                        print(
                            "{:.2f} GB remaining in the output directory {}."
                            .format(free_space, self.output_directory)
                        )
                        while True:
                            answer = input(
                                "Do you want to continue? [Y/n] "
                            ).lower()
                            if answer == "y" or answer == "":
                                print("Ok, disabling further disk space checks.")
                                self.disk_limit = 0
                                break
                            elif answer == "n":
                                return True  # stop
                            else:
                                print('Unknown answer.  Please answer with "y" or "n".')
        return False  # no problem

    def finalize(self, state, t, n):
        """Write the final state to the history file if necessary.

        Arguments:
         - state: final state of the model run
         - t: last point in time of the model run
         - n: total number of integration steps in the model run
        """
        # Write data to the history file if it is time to do so
        if n != self.last_saved_frame:
            self.write_history_file(state, t, n)

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

    def create_history_file(self, grid: Grid):
        """Create a new history file."""
        with nc.Dataset(self.hist_path, "w") as ncfile:
            # Store the experiment parameters
            ncfile.setncatts(self.experiment_parameters)

            # Create the dimensions
            ncfile.createDimension("t")  # unlimited size
            ncfile.createDimension("x_b", grid.nx)
            ncfile.createDimension("y_b", grid.ny)
            ncfile.createDimension("z_b", grid.nz)
            ncfile.createDimension("x_u", grid.nx)
            ncfile.createDimension("y_u", grid.ny)
            ncfile.createDimension("z_u", grid.nz)
            ncfile.createDimension("x_v", grid.nx)
            ncfile.createDimension("y_v", grid.ny)
            ncfile.createDimension("z_v", grid.nz)
            ncfile.createDimension("x_w", grid.nx)
            ncfile.createDimension("y_w", grid.ny)
            ncfile.createDimension("z_w", grid.nz)
            ncfile.createDimension("x_vor_i", grid.nx)
            ncfile.createDimension("y_vor_i", grid.ny)
            ncfile.createDimension("z_vor_i", grid.nz)
            ncfile.createDimension("x_vor_j", grid.nx)
            ncfile.createDimension("y_vor_j", grid.ny)
            ncfile.createDimension("z_vor_j", grid.nz)
            ncfile.createDimension("x_vor_k", grid.nx)
            ncfile.createDimension("y_vor_k", grid.ny)
            ncfile.createDimension("z_vor_k", grid.nz)

            # Create the variables with one dimension
            v = ncfile.createVariable("n", int, ("t",))
            v.long_name = "integration step in the model run"

            v = ncfile.createVariable("t", float, ("t",))
            v.long_name = "time in the model run"
            v.units = self.unit("T")

            v = ncfile.createVariable("x_b", float, ("x_b",))
            v.long_name = grid.x_b.name
            v.units = self.unit(grid.x_b.dimension)
            v[:] = grid.x_b_1D
            v = ncfile.createVariable("y_b", float, ("y_b",))
            v.long_name = grid.y_b.name
            v.units = self.unit(grid.y_b.dimension)
            v[:] = grid.y_b_1D
            v = ncfile.createVariable("z_b", float, ("z_b",))
            v.long_name = grid.z_b.name
            v.units = self.unit(grid.z_b.dimension)
            v[:] = grid.z_b_1D

            v = ncfile.createVariable("x_u", float, ("x_u",))
            v.long_name = grid.x_vel["i"].name
            v.units = self.unit(grid.x_vel["i"].dimension)
            v[:] = grid.x_u_1D
            v = ncfile.createVariable("y_u", float, ("y_u",))
            v.long_name = grid.y_vel["i"].name
            v.units = self.unit(grid.y_vel["i"].dimension)
            v[:] = grid.y_u_1D
            v = ncfile.createVariable("z_u", float, ("z_u",))
            v.long_name = grid.z_vel["i"].name
            v.units = self.unit(grid.z_vel["i"].dimension)
            v[:] = grid.z_u_1D

            v = ncfile.createVariable("x_v", float, ("x_v",))
            v.long_name = grid.x_vel["j"].name
            v.units = self.unit(grid.x_vel["j"].dimension)
            v[:] = grid.x_v_1D
            v = ncfile.createVariable("y_v", float, ("y_v",))
            v.long_name = grid.y_vel["j"].name
            v.units = self.unit(grid.y_vel["j"].dimension)
            v[:] = grid.y_v_1D
            v = ncfile.createVariable("z_v", float, ("z_v",))
            v.long_name = grid.z_vel["j"].name
            v.units = self.unit(grid.z_vel["j"].dimension)
            v[:] = grid.z_v_1D

            v = ncfile.createVariable("x_w", float, ("x_w",))
            v.long_name = grid.x_vel["k"].name
            v.units = self.unit(grid.x_vel["k"].dimension)
            v[:] = grid.x_w_1D
            v = ncfile.createVariable("y_w", float, ("y_w",))
            v.long_name = grid.y_vel["k"].name
            v.units = self.unit(grid.y_vel["k"].dimension)
            v[:] = grid.y_w_1D
            v = ncfile.createVariable("z_w", float, ("z_w",))
            v.long_name = grid.z_vel["k"].name
            v.units = self.unit(grid.z_vel["k"].dimension)
            v[:] = grid.z_w_1D

            v = ncfile.createVariable("x_vor_i", float, ("x_vor_i",))
            v.long_name = grid.x_vor["i"].name
            v.units = self.unit(grid.x_vor["i"].dimension)
            v[:] = grid.x_vor_i_1D
            v = ncfile.createVariable("y_vor_i", float, ("y_vor_i",))
            v.long_name = grid.y_vor["i"].name
            v.units = self.unit(grid.y_vor["i"].dimension)
            v[:] = grid.y_vor_i_1D
            v = ncfile.createVariable("z_vor_i", float, ("z_vor_i",))
            v.long_name = grid.z_vor["i"].name
            v.units = self.unit(grid.z_vor["i"].dimension)
            v[:] = grid.z_vor_i_1D

            v = ncfile.createVariable("x_vor_j", float, ("x_vor_j",))
            v.long_name = grid.x_vor["j"].name
            v.units = self.unit(grid.x_vor["j"].dimension)
            v[:] = grid.x_vor_j_1D
            v = ncfile.createVariable("y_vor_j", float, ("y_vor_j",))
            v.long_name = grid.y_vor["j"].name
            v.units = self.unit(grid.y_vor["j"].dimension)
            v[:] = grid.y_vor_j_1D
            v = ncfile.createVariable("z_vor_j", float, ("z_vor_j",))
            v.long_name = grid.z_vor["j"].name
            v.units = self.unit(grid.z_vor["j"].dimension)
            v[:] = grid.z_vor_j_1D

            v = ncfile.createVariable("x_vor_k", float, ("x_vor_k",))
            v.long_name = grid.x_vor["k"].name
            v.units = self.unit(grid.x_vor["k"].dimension)
            v[:] = grid.x_vor_k_1D
            v = ncfile.createVariable("y_vor_k", float, ("y_vor_k",))
            v.long_name = grid.y_vor["k"].name
            v.units = self.unit(grid.y_vor["k"].dimension)
            v[:] = grid.y_vor_k_1D
            v = ncfile.createVariable("z_vor_k", float, ("z_vor_k",))
            v.long_name = grid.z_vor["k"].name
            v.units = self.unit(grid.z_vor["k"].dimension)
            v[:] = grid.z_vor_k_1D

            # TODO: add mask if a mask is implemented

            # Create variables for the model data.
            for variable in self.hist_variables:
                if variable.get_nature() == "scalar":
                    v = ncfile.createVariable(
                        variable.nickname, float, ("t", "z_b", "y_b", "x_b")
                    )
                    v.long_name = variable.name
                    v.units = self.unit(variable.dimension)
                elif variable.get_nature() == "velocity":
                    for i, u in zip("ijk", "uvw"):
                        v = ncfile.createVariable(
                            variable[i].nickname, float,
                            ("t", "z_" + u, "y_" + u, "x_" + u)
                        )
                        v.long_name = variable[i].name
                        v.units = self.unit(variable[i].dimension)
                elif variable.get_nature() == "vorticity":
                    for i in "ijk":
                        v = ncfile.createVariable(
                            variable[i].nickname, float,
                            ("t", "z_vor_" + i, "y_vor_" + i, "x_vor_" + i)
                        )
                        v.long_name = variable[i].name
                        v.units = self.unit(variable[i].dimension)

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
        self.last_saved_frame = n

    def get_disk_space_in_GB(self):
        """Return the available disk space of the output directory in GB.

        Explanation: https://stackoverflow.com/a/12327880/3661532
        """
        statvfs = os.statvfs(self.output_directory)
        return statvfs.f_frsize * statvfs.f_bavail / 1e9


if __name__ == "__main__":
    import numpy as np

    from variables import get_state


    param = {
        ### Parameters necessary for the IO class
        "datadir": "~/data/Nyles",
        "expname": "test_exp",
        "timestep_history": 1.0,
        "disk_space_warning": 0.5,
        "unit_length": "m",
        "unit_duration": "s",
        ## Choose between a list or "all"
        "variables_in_history": "all",
        # "variables_in_history": ["u", "b"],
        ## Select one of the following options
        "mode": "overwrite",
        # "mode": "count",
        # "mode": "continue",
        ### Parameters necessary for State and/or Grid
        "Lx": 60,
        "Ly": 50,
        "Lz": 10,
        "nx": 6,
        "ny": 5,
        "nz": 4,
        "nh": 3,
        "neighbours": {},  # TODO: try with neighbours when halo-handling is implemented
        ### Parameters to test the storage of parameters in the history file
        "a boolean variable": True,
        "a 2nd boolean variable": False,
        "a long list": list(range(50)),
    }
    state = get_state(param)
    print("* Full state:")
    print(state)

    grid = Grid(param)

    io = NylesIO(param)
    print("* Output directory:", io.output_directory)
    print("* History file:", io.hist_path)
    print("* Experiment parameters to save:")
    print(io.experiment_parameters)
    io.init(state, grid)
    print("* Variables in the history file:")
    print(io.hist_variables)
    b = state.get("b").view("i")
    for i, t in enumerate([0.2, 0.5, 0.7, 1.0, 1.5, 2.1, 3.0]):
        n = i+1
        print("* t = {:9.2f} / n = {:6} / n_hist = {:6}".format(t, n, io.n_hist))
        b[:] = t * np.ones_like(b)
        io.do(state, t, n)
