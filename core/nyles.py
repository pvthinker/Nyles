"""The Nyles main class."""

import numpy as np

import model_les_AL as model_LES
import model_advection as model_adv
import variables
import grid
import nylesIO


class Nyles(object) :
    """
    Attributes :
        - self.grid
        - self.model
        - self.IO
        - self.tend
        - self.auto_dt: boolean
        - self.cfl
        - self.dt0

    Methods :
        *private :
            - initiate(param) : loads the desired model and initiates the IO
            defines tend and cfl from the param
            - compute_dt(): calculate the timestep
        *public :
            - run() : main loop. Iterates the model forward and saves the state
            in the history file
    """
    def __init__(self, param):
        self.grid = grid.Grid(param) #loads the grid
        self.IO = nylesIO.NylesIO(param) #loads the IO
        self.initiate(param) #initiates the model and needed variables

    def initiate(self, param):
        if param['modelname'] == 'LES' :
            self.model = model_LES.LES(param)
        elif param['modelname'] == 'adv' :
            self.model = model_adv.Advection(param)

        self.tend = param['tend']
        self.auto_dt = bool(param['auto_dt'])  # enforce boolean type
        self.cfl = param['cfl']
        self.dt0 = param['dt']

    def run(self):
        t = 0.0
        n = 0

        print("Creating output file:", self.IO.hist_path)
        self.IO.init(self.model.state, self.grid, t, n)

        while True:
            dt = self.compute_dt()
            self.model.forward(t, dt)
            t += dt
            n += 1
            self.IO.do(self.model.state, t, n)
            if t >= self.tend:
                break

        self.IO.finalize(self.model.state, t, n)

    def compute_dt(self):
        """Calculate timestep dt from contravariant velocity U and cfl.

        The fixed value self.dt0 is returned if and only if self.auto_dt
        is False or the velocity is zero everywhere (this is an
        extremely rare case, but must be included to avoid division by
        zero).  Otherwise, the timestep is calculated as
            dt = cfl / max(|U|) .
        Note that the "dx" is hidden in the contravariant velocity U,
        which has units of 1/T.  In the formula, |U| denotes the norm of
        the contravariant velocity vector.
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
                return self.dt0
            else:
                return self.cfl / U_max
        else:
            return self.dt0


if __name__ == "__main__":
    import topology as topo

    procs = [1, 1, 1]
    topo.topology = "closed"
    myrank = 0
    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    param = {
        # IO parameters
        "datadir": "~/data/Nyles",
        "expname": "test_exp",
        "timestep_history": 1.0,
        "variables_in_history": "all",
        "mode": "overwrite",
        # Grid parameters
        "Lx": 1.0,
        "Ly": 1.0,
        "Lz": 1.0,
        "nx": 32,
        "ny": 32,
        "nz": 16,
        "nh": 0,
        # Other parameters
        "neighbours": neighbours,
        "modelname": "LES",
        "timestepping": "LFAM3",
        "tend": 10.0,
        "auto_dt": True,
        "cfl": 1.0,
        "dt": 0.2,
    }

    nyles = Nyles(param)
    U = nyles.model.state.U["i"].view()
    V = nyles.model.state.U["j"].view()
    W = nyles.model.state.U["k"].view()
    if param["auto_dt"]:
        print("Case 1: U = 0, V = 0, W = 0")
        print("    dt is", nyles.compute_dt())
        print("should be", param["dt"])
        U += 1
        print("Case 2: U = 1, V = 0, W = 0")
        print("    dt is", nyles.compute_dt())
        print("should be", param["cfl"] / np.sqrt(1))
        V += 1
        print("Case 3: U = 1, V = 1, W = 0")
        print("    dt is", nyles.compute_dt())
        print("should be", param["cfl"] / np.sqrt(2))
        W += 1
        print("Case 4: U = 1, V = 1, W = 1")
        print("    dt is", nyles.compute_dt())
        print("should be", param["cfl"] / np.sqrt(3))
