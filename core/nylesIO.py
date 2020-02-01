"""Provide Input/Output handling for Nyles through the class NylesIO.

The NylesIO class serves for two closely related purposes:
 1. Output: Store the output of a model run in a netCDF file.
 2. Input:  Load the output of a previously saved model run and continue
            the experiment from the last saved point.
It also provides methods that are useful for the handling of other
output files.

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

# Standard library imports
import os
import shutil
from enum import Enum   # since Python 3.4

# Third party import
import netCDF4 as nc

# Local imports
from variables import State
from grid import Grid
import mpitools


class NylesIO(object):
    """Input/Output handler for Nyles.

    The interaction with an instance of this class should usually only
    occur via the public attributes and methods listed below, and not by
    accessing the private attributes or methods.

    Attributes for public access:
     - output_directory: path to the output directory
     - hist_path: path to the history file; this string contains the
        path to the output directory
     - script_path: path of the script file backup

    Attributes mostly for private access:
     - variables_in_history: list or string describing the variables
        that the user wants to save in the history file
     - hist_variables: dictionary containing all the scalars that
        are actually saved in the history file; this is created on
        “init” together with the history file; the keys are the names
        used in the history file, the values are the nicknames used in
        the model
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
     - backup_scriptfile

    Methods mostly for private access:
     - create_history_file
     - write_history_file
     - get_disk_space_in_GB
    """

    MAX_LENGTH_ATTRIBUTE = 100

    def __init__(self, param):
        """Create an object to handle Input and Output.

        The argument “param” determines the behaviour of this object.
        Furthermore, every key-value-pair in "param" is saved in the
        history file.  Values of type int, float or string are saved as
        they are.  Boolean variables are converted to the string "True"
        or "False".  Every other variable is converted to a string and
        shortened to MAX_LENGTH_ATTRIBUTE characters if necessary.
        """
        self.disk_limit = param["disk_space_warning"]
        self.simplified_grid = param["simplified_grid"]
        self.variables_in_history = param["variables_in_history"]
        self.dt_hist = param["timestep_history"]
        self.t_next_hist = 0.0
        self.n_hist = 0
        self.last_saved_frame = None
        self.myrank = param["myrank"]
        # Define a function to replace dimensions by units for convenience
        self.unit = lambda dimension: (
            dimension.replace("T", param["unit_duration"]).replace("L", param["unit_length"])
        )
        self.n_tracers = param["n_tracers"]

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
                    "{}: {}".format(type(value), repr(value))
                    [:self.MAX_LENGTH_ATTRIBUTE+1]
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
        self.hist_path = os.path.join(out_dir, expname +
                                      "_%02i" % param["myrank"] +
                                      "_hist.nc")
        self.script_path = os.path.join(out_dir, expname + ".py")
        self.output_directory = out_dir

    def init(self, state: State, grid: Grid, t=0.0, n=0):
        """Create a new history file with the initial data.

        Arguments:
         - state: initial state of the model run
         - t: initial time of the model run, usually 0.0
         - n: integration step of the model, usually 0
        """
        # Check the list of variables to save in the history file
        if self.variables_in_history == "all":
            self.variables_in_history = state.toc.keys()
        elif self.variables_in_history == "prognostic":
            self.variables_in_history = state.get_prognostic_variables()
        elif self.variables_in_history == "p+p":
            self.variables_in_history = state.get_prognostic_variables()
            if "p" in state.toc:
                self.variables_in_history.append("p")
            else:
                raise ValueError(
                    "pressure is not a variable of this model; "
                    "please change the parameter variables_in_history."
                )
        elif any(variable not in state.toc for variable in self.variables_in_history):
            raise ValueError(
                "a variable chosen for the history file does not exist in this "
                "model.  The available variables are " + str(state.toc)
            )
        else:
            # Add tracers
            for i in range(self.n_tracers):
                nickname = "t{}".format(i)
                if nickname not in self.variables_in_history:
                    self.variables_in_history.append(nickname)

        # Check that the user did not forget to add variables
        if len(self.variables_in_history) == 0:
            raise ValueError(
                "no variables are selected for the history file.  "
                "This means, all calculation results are lost.  "
                "I refuse to waste energy like this."
            )

        # Make list of strings into list of Scalar and Vector objects
        variables = [state.get(v) for v in self.variables_in_history]

        # Create the output directory if necessary
        if not os.path.isdir(self.output_directory):
            if self.myrank == 0:
                os.makedirs(self.output_directory)
        # block all ranks until rank 0 has created the folder(s)
        mpitools.barrier()
        # let be very cautious
        assert os.path.isdir(self.output_directory)

        # Create the history file and save the initial state
        self.create_history_file(grid, variables)
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
        if t < self.t_next_hist:
            # Nothing to do yet
            return False
        # Otherwise it is time to write data to the history file
        self.write_history_file(state, t, n)
        self.t_next_hist += self.dt_hist
        if self.disk_limit <= 0:
            # No disk space limit defined
            return False
        # Otherwise check if the remaining disk space is sufficient
        try:
            free_space = self.get_disk_space_in_GB()
        except Exception as e:
            print("")
            print("Cannot determine available disk space.")
            print("Error message:", e)
            print("Disabling further disk space checks.")
            self.disk_limit = 0
            # This is not a reason to stop the program
            return False
        if free_space >= self.disk_limit:
            # Everything fine
            return False
        # Otherwise, print a warning
        print("")
        print("!"*50)
        print("Warning, low disk space:")
        print(
            "{:.2f} GB remaining in the output directory {}"
            .format(free_space, self.output_directory)
        )
        # Ask the user what to do
        while True:
            answer = input(
                "Do you want to continue? [Y/n] "
            ).lower()
            if answer == "n":
                # Stop the program
                return True
            elif answer == "y" or answer == "":
                # Check if more space is free now
                try:
                    free_space = self.get_disk_space_in_GB()
                except:
                    # Unknown error; check free space again next time
                    pass
                else:
                    # elif cannot be used here
                    if free_space < self.disk_limit:
                        # The problem persists
                        print("Disabling further disk space checks.")
                        self.disk_limit = 0
                # Continue the program
                return False
            else:
                print('Unknown answer.', end='  ')
                print('Please answer with "y" or "n".')

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

    def create_history_file(self, grid: Grid, variables: list):
        """Create a new history file with “variables” living on “grid”."""
        # Open new netCDF file for writing data to it (or for
        # overwriting existing data)
        with nc.Dataset(self.hist_path, "w") as ncfile:
            # Store the experiment parameters
            ncfile.setncatts(self.experiment_parameters)

            # Create the dimensions
            ncfile.createDimension("t")  # unlimited size
            for x, i in zip("xyz", "ijk"):
                if self.simplified_grid:
                    ncfile.createDimension("{}".format(x), grid.size[i])
                    continue
                # “p” stands for “point”: b-point, u-point, v-point, etc.
                for p in ["b", "u", "v", "w", "vor_i", "vor_j", "vor_k"]:
                    ncfile.createDimension("{}_{}".format(x, p), grid.size[i])

            # Create the variables with one dimension
            v = ncfile.createVariable("n", int, ("t",))
            v.long_name = "integration step in the model run"

            v = ncfile.createVariable("t", float, ("t",))
            v.long_name = "time in the model run"
            v.units = self.unit("T")

            if self.simplified_grid:
                v = ncfile.createVariable("x", float, ("x",))
                v.long_name = grid.x_b.name
                v.units = self.unit(grid.x_b.dimension)
                v[:] = grid.x_b_1D
                v = ncfile.createVariable("y", float, ("y",))
                v.long_name = grid.y_b.name
                v.units = self.unit(grid.y_b.dimension)
                v[:] = grid.y_b_1D
                v = ncfile.createVariable("z", float, ("z",))
                v.long_name = grid.z_b.name
                v.units = self.unit(grid.z_b.dimension)
                v[:] = grid.z_b_1D

            else:
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

            # Create variables for the model data and create a dict of
            # the history file variables to make writing new data easy
            self.hist_variables = {}
            for variable in variables:
                nickname = variable.nickname
                nature = variable.get_nature()
                if nature == "scalar":
                    # Use in the history file the same name as in the model
                    hist_name = nickname
                    v = ncfile.createVariable(
                        hist_name,
                        float,
                        ("t", "z", "y", "x") if self.simplified_grid else
                        ("t", "z_b", "y_b", "x_b"),
                    )
                    v.long_name = variable.name
                    v.units = self.unit(variable.dimension)
                    self.hist_variables[hist_name] = nickname
                elif nature == "velocity":
                    if nickname == "u":
                        # Use in the history file the name u/v/w
                        modifier = str.lower
                    elif nickname == "U":
                        # Use in the history file the name U/V/W
                        modifier = str.upper
                    else:
                        raise NotImplementedError(
                            "unknown kind of velocity: " + nickname
                        )
                    for i, u in zip("ijk", "uvw"):
                        hist_name = modifier(u)
                        v = ncfile.createVariable(
                            hist_name,
                            float,
                            ("t", "z", "y", "x") if self.simplified_grid else
                            ("t", "z_" + u, "y_" + u, "x_" + u),
                        )
                        v.long_name = variable[i].name
                        v.units = self.unit(variable[i].dimension)
                        self.hist_variables[hist_name] = variable[i].nickname
                elif nature == "vorticity":
                    for i in "ijk":
                        # Use in the history file the same name as in the model
                        hist_name = variable[i].nickname
                        v = ncfile.createVariable(
                            hist_name,
                            float,
                            ("t", "z", "y", "x") if self.simplified_grid else
                            ("t", "z_vor_" + i, "y_vor_" + i, "x_vor_" + i),
                        )
                        v.long_name = variable[i].name
                        v.units = self.unit(variable[i].dimension)
                        self.hist_variables[hist_name] = variable[i].nickname
                else:
                    raise ValueError(
                        "unknown nature", nature, "of variable", nickname
                    )

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
            for hist_name, nickname in self.hist_variables.items():
                ncfile[hist_name][self.n_hist] = state.get(nickname).view("i")
        self.n_hist += 1
        self.last_saved_frame = n

    def get_disk_space_in_GB(self):
        """Return the available disk space of the output directory in GB.

        Explanation: https://stackoverflow.com/a/12327880/3661532
        """
        statvfs = os.statvfs(self.output_directory)
        return statvfs.f_frsize * statvfs.f_bavail / 1e9

    def backup_scriptfile(self, filename):
        shutil.copyfile(filename, self.script_path)

    def write_githashnumber(self):
        gitfile = self.output_directory +'/nyles.git'
        import version
        githash = version.get_nyles_hash_number()
        # print(gitfile)
        with open(gitfile, 'w') as fid:
            fid.write('# this experiment has been done with\n')
            fid.write('# Nyles commit\n')
            fid.write(githash+'\n')
            fid.write('# to rerun it with same version \n')
            fid.write('# git checkout %s' % githash[:7])

if __name__ == "__main__":
    import numpy as np

    from variables import get_state

    param = {
        # Parameters necessary for the IO class
        "datadir": "~/data/Nyles",
        "expname": "test_exp",
        "timestep_history": 1.0,
        "disk_space_warning": 0.5,
        "unit_length": "m",
        "unit_duration": "s",
        "n_tracers": 0,  # does not change anything here
        "simplified_grid": False,
        "variables_in_history": "p+p",
        # Select one of the following options
        "mode": "overwrite",
        # "mode": "count",
        # "mode": "continue",
        # Parameters necessary for State and/or Grid
        "Lx": 60,
        "Ly": 50,
        "Lz": 10,
        "nx": 6,
        "ny": 5,
        "nz": 4,
        "nh": 3,
        "neighbours": {},  # TODO: try with neighbours when halo-handling is implemented
        # Parameters to test the storage of parameters in the history file
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
        b[:] = t * np.ones_like(b)
        io.do(state, t, n)
        print("* t = {:9.2f} / n = {:6} / n_hist = {:6}".format(t, n, io.n_hist))
