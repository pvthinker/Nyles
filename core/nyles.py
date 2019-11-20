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
            - compute_dt()
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
        """Calculate the timestep as function of velocity and grid size.

        This function returns the fixed value self.dt0 if self.auto_dt
        is False or if the velocity is everywhere zero.  Otherwise it
        calculates the current timestep dt using self.cfl, the size of
        the grid and the maximal velocity.  The formula (in x) is:
            dt = cfl * dx^2 / u
        Note that dx is squared, because u is the covariant velocity,
        which has dimensions of L^2/T."""
        if self.auto_dt:
            # Compute maximal velocities
            u_max = np.max(np.abs(self.state.u["i"].view()))
            v_max = np.max(np.abs(self.state.u["j"].view()))
            w_max = np.max(np.abs(self.state.u["k"].view()))
            # If all velocities are zero, take the standard timestep
            if u_max == 0.0 and v_max == 0.0 and w_max == 0.0:
                return self.dt0
            # Otherwise compute dt with cfl, the grid size, and the velocities
            dt_x = self.cfl * self.grid.dx2 / u_max
            dt_y = self.cfl * self.grid.dy2 / v_max
            dt_z = self.cfl * self.grid.dz2 / w_max
            # Return the smallest non-zero value of dt
            return np.nanmin([
                dt_x if dt_x else np.nan,
                dt_y if dt_y else np.nan,
                dt_z if dt_z else np.nan,
            ])
        else:
            return self.dt0
