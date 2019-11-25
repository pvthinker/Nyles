"""The Nyles main class."""

import numpy as np

import model_les_AL as model_LES
import model_advection as model_adv
import variables
import grid
import nylesIO
import plotting
from timing import timing


class Nyles(object) :
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
        # Check and get user parameters
        user_param.check()
        param = user_param.view_parameters()

        # Load the IO; only the parameters modifiable by the user are saved
        self.IO = nylesIO.NylesIO(param)

        # Add parameters that are automatically set
        param["nx"] = param["global_nx"] // param["npx"]
        param["ny"] = param["global_ny"] // param["npy"]
        param["nz"] = param["global_nz"] // param["npz"]
        # TODO: make the following general
        # for reference (taken from an experiment file; remove when done):
        # import topology as topo
        # topo.topology = 'closed'
        # procs = [1, 1, 1]
        # myrank = 0
        # loc = topo.rank2loc(myrank, procs)
        # neighbours = topo.get_neighbours(loc, procs)
        param.update({
            "procs": [1, 1, 1], "neighbours": {},
            "topology": "closed",  # is this always the same as param["geometry"]?
            "npre": 3, "npost": 3, "omega": 0.8, "ndeepest": 20, "maxite": 20,
            "tol": 1e-6,
        })

        # Load the grid with the extended parameters
        self.grid = grid.Grid(param)

        # Initiate the model and needed variables
        self.initiate(param)

    def initiate(self, param):
        if param['modelname'] == 'LES' :
            self.model = model_LES.LES(param, self.grid)
        elif param['modelname'] == 'advection':
            self.model = model_adv.Advection(param, self.grid)

        self.tend = param['tend']
        self.auto_dt = param['auto_dt']
        self.dt0 = param['dt']
        self.cfl = param['cfl']
        self.dt_max = param['dt_max']

        # Load the plotting module
        self.plotting = plotting.Plotting(param, self.model.state, self.grid)

    def run(self):
        t = 0.0
        n = 0

        print("Creating output file:", self.IO.hist_path)
        self.IO.init(self.model.state, self.grid, t, n)

        # Open the plotting window and draw the initial state
        self.plotting.init(t, n)
        # Give the user some time to enjoy the first frame.
        # The rest of the simulation won't get prettier,
        # as long as we haven't fixed the numerics.
        print("Resize the window to a suitable size,")
        print("move the camera into a good angle,")
        print("lean back in your seat and ...")
        input("... press Enter to start! ")

        while True:
            dt = self.compute_dt()
            self.model.forward(t, dt)
            t += dt
            n += 1
            self.IO.do(self.model.state, t, n)
            self.plotting.update(t, n)
            if t >= self.tend:
                break

        self.IO.finalize(self.model.state, t, n)

    @timing
    def compute_dt(self):
        """Calculate timestep dt from contravariant velocity U and cfl.

        The fixed value self.dt0 is returned if and only if self.auto_dt
        is False.  Otherwise, the timestep is calculated as
            dt = cfl / max(|U|) .
        Note that the "dx" is hidden in the contravariant velocity U,
        which has units of 1/T.  In the formula, |U| denotes the norm of
        the contravariant velocity vector.  If the velocity is zero
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
            #
            # One can speed up the calculation by using instead
            # np.max([np.max(np.abs(U)), np.max(np.abs(V)), np.max(np.abs(W))])
            # to approximate U_max.  However, this gives only a lower
            # bound that can differ from the actual value by a factor of
            # up to sqrt(3) = 1.73.  So with this approximation, the
            # timestep dt is chosen too big, which could lead to an
            # unstable integration if cfl is not adjusted accordingly.
            U_max = np.sqrt(np.max(U**2 + V**2 + W**2))
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
