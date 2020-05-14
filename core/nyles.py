"""The Nyles main class."""

import sys

import numpy as np
import pickle

import model_les
import model_advection as model_adv
import variables
import grid
import nylesIO
import plotting
import timing
import topology as topo
import mpitools
from time import time


class Nyles(object):
    """
    Attributes :
        - self.grid
        - self.model
        - self.IO
        - self.tend
        - self.auto_dt
        - self.dt0
        - self.cfl
        - self.dt_max

    Methods :
        *private :
            - initiate(param): loads the desired model and copies the
            time variables from param
            - compute_dt(): calculate the timestep
        *public :
            - run() : main loop. Iterates the model forward and saves the state
            in the history file
    """

    def __init__(self, user_param):
        # Check, freeze, and get user parameters
        user_param.check()
        user_param.freeze()
        param = user_param.view_parameters()
        self.param = param
        npx = param["npx"]
        npy = param["npy"]
        npz = param["npz"]

        # Add parameters that are automatically set
        param["nx"] = param["global_nx"] // npx
        param["ny"] = param["global_ny"] // npy
        param["nz"] = param["global_nz"] // npz

        # Set up MPI
        topo.topology = param["geometry"]
        procs = [npz, npy, npx]
        myrank = mpitools.get_myrank(procs)
        loc = topo.rank2loc(myrank, procs)
        neighbours = topo.get_neighbours(loc, procs)
        param["procs"] = procs
        param["myrank"] = myrank
        param["neighbours"] = neighbours
        param["loc"] = loc
        self.myrank = myrank

        if self.myrank == 0:
            print("-"*80)
        self.banner()

        # Load the grid with the extended parameters
        self.grid = grid.Grid(param)

        # Load the IO; only the parameters modifiable by the user are saved
        self.IO = nylesIO.NylesIO(param)

        # redirect the x11 to output.txt
        #sys.stdout = Logger(self.IO.output_directory+'/output.txt')

        # save the param dictionary in a pkl file
        if self.myrank == 0:
            with open(self.IO.output_directory+"/param.pkl", "wb") as fid:
                pickle.dump(param, fid)

        # Initiate the model and needed variables
        self.initiate(param)

    def initiate(self, param):
        if param['modelname'] == 'LES':
            self.model = model_les.LES(param, self.grid)
        elif param['modelname'] == 'linear':
            self.model = model_les.LES(param, self.grid, linear=True)
        elif param['modelname'] == 'advection':
            self.model = model_adv.Advection(param, self.grid)

        self.tend = param['tend']
        self.auto_dt = param['auto_dt']
        self.dt0 = param['dt']
        self.cfl = param['cfl']
        self.dt_max = param['dt_max']

        # Load the plotting module
        if param["show"]:
            self.plotting = plotting.Plotting(
                param, self.model.state, self.grid)
        else:
            self.plotting = None
        # number of grid cell per subdomain
        self.gridcellpersubdom = param["nx"]*param["ny"]*param["nz"]

    def run(self):
        t = 0.0
        n = 0
        self.model.diagnose_var(self.model.state)

        # Open the plotting window and draw the initial state
        if self.plotting:
            self.plotting.init(t, n)
            print("Resize the window to a suitable size,")
            print("move the camera into a good angle,")
            print("lean back in your seat and ...")
            input("... press Enter to start! ")

        if self.myrank == 0:
            print("Creating output file:", self.IO.hist_path)

        self.IO.init(self.model.state, self.grid, t, n)

        if self.myrank == 0:
            print("Backing up script to:", self.IO.script_path)
            self.IO.backup_scriptfile(sys.argv[0])
            try:
                self.IO.write_githashnumber()
            except:
                # this occurs at TGCC where the code
                # can't be installed from github
                print("Failed to save the git version of the code experiment")

        time_length = len(str(int(self.tend))) + 3
        time_string = "\r"+", ".join([
            "n = {:3d}",
            "t = {:" + str(time_length) + ".2f}/{:" +
            str(time_length) + ".2f}",
            "dt = {:.4f}",
            "perf = {:.2e}"
        ])

        if self.myrank == 0:
            print("-"*80)
        stop = False

        realtime0 = time()

        while not(stop):
            dt = self.compute_dt()
            blowup = self.model.forward(t, dt)
            t += dt
            n += 1
            stop = self.IO.write(self.model.state, t, n)
            if self.myrank == 0:
                realtime = time()
                # cpu per iteration per grid cell (in second)
                # this is very good metric of how fast the code runs
                # typical values are O(4e-6)
                # > 1e-5 is not good (e.g. too large subdomains)
                # < 2e-6 is excellent, report it to the developpers!
                perf = (realtime-realtime0)/(self.gridcellpersubdom)
                realtime0 = realtime
                print(time_string.format(n, t, self.tend, dt, perf), end='')
            # Blowup might be detected on one subdomain
            # and not on the other.
            # We need to spread the word that blowup
            # had occured and that all cores should stop.
            # An easy way to do it is to do a global sum
            # of localmood (0=fine, 1=blowup)
            localmood = 1. if blowup else 0.
            globalmood = mpitools.global_sum(localmood)
            # and test globalmood
            if globalmood > 0.:
                if self.myrank == 0:
                    print('')
                    print('BLOW UP!')
                    print('write a last snapshot at blowup time', end="")
                # force a last write in the netCDF
                # might help locate the origin of the
                # blowup
                self.IO.t_next_hist = t
                stop = self.IO.write(self.model.state, t, n)

            if self.plotting:
                self.plotting.update(t, n)
            if t >= self.tend or stop:
                break
            mpitools.barrier()
        if self.myrank == 0:
            print()
            print("-"*80)
        if self.myrank == 0:
            if stop:
                print("Job is aborted")
            else:
                print("Job completed as expected")

        self.IO.finalize(self.model.state, t, n)

        if self.myrank == 0:
            print("Output written to:", self.IO.hist_path)

        self.model.write_stats(self.IO.output_directory)

        timing.write_timings(self.IO.output_directory)
        timing.analyze_timing(self.IO.output_directory)

    def compute_dt(self):
        """Calculate timestep dt from contravariant velocity U and cfl.

        The fixed value self.dt0 is returned if and only if self.auto_dt
        is False.  Otherwise, the timestep is calculated as
            dt = cfl / max(|U|) .
        Note that the "dx" is hidden in the contravariant velocity U,
        which has dimension 1/T.  In the formula, |U| denotes the norm
        of the contravariant velocity vector.  If the velocity is zero
        everywhere, self.dt_max is returned.  Otherwise, the smaller
        value of dt and dt_max is returned.
        """
        if self.auto_dt:
            # Get U, V, and W in the same orientation
            U_object = self.model.state.U["i"]
            U = U_object.view()
            V = self.model.state.U["j"].viewlike(U_object)
            W = self.model.state.U["k"].viewlike(U_object)
            # Since the sqrt-function is strictly monotonically
            # increasing, the order of sqrt and max can be exchanged.
            # This way, it is only necessary to calculate the square
            # root of a single value, which is faster.
            U_max = np.sqrt(np.max(U**2 + V**2 + W**2))
            U_max = mpitools.global_max(U_max)
            # Note: the if-statement cannot be replaced by try-except,
            # because U_max is a numpy-float which throws a warning
            # instead of an error in case of division by zero.
            if U_max == 0.0:
                return self.dt_max
            else:
                dt = self.cfl / U_max
                return min(dt, self.dt_max)
        else:
            return self.dt0

    def banner(self):
        logo = [
            "  _   _       _            ",
            " | \ | |     | |           ",
            " |  \| |_   _| | ___  ___  ",
            " | . ` | | | | |/ _ \/ __| ",
            " | |\  | |_| | |  __/\__ \ ",
            " |_| \_|\__, |_|\___||___/ ",
            "         __/ |             ",
            "        |___/              ",
            "                           "]
        if self.myrank == 0:
            print("Welcome to")
            for l in logo:
                print(" "*10+l)


class Logger(object):
    def __init__(self, logfile):
        self.terminal = sys.stdout
        self.log = open(logfile, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


if __name__ == "__main__":
    from parameters import UserParameters

    param = UserParameters()

    param.time["dt_max"] = 1.5

    nyles = Nyles(param)
    U = nyles.model.state.U["i"].view()
    V = nyles.model.state.U["j"].view()
    W = nyles.model.state.U["k"].view()
    if param.time["auto_dt"]:
        print("Case 1: U = 0, V = 0, W = 0")
        print("    dt is", nyles.compute_dt())
        print("should be", param.time["dt_max"])
        U += 1
        print("Case 2: U = 1, V = 0, W = 0")
        print("    dt is", nyles.compute_dt())
        print("should be", param.time["cfl"] / np.sqrt(1))
        V += 1
        print("Case 3: U = 1, V = 1, W = 0")
        print("    dt is", nyles.compute_dt())
        print("should be", param.time["cfl"] / np.sqrt(2))
        W += 1
        print("Case 4: U = 1, V = 1, W = 1")
        print("    dt is", nyles.compute_dt())
        print("should be", param.time["cfl"] / np.sqrt(3))
