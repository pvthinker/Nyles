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
        - self.cfl
        - self.dx
        - self.dt
    dx will be an attribute of the grid

    Methods :
        *private :
            - __initiate(param) : loads the desired model and initiates the IO
            defines tend and cfl from the param
        *public :
            - run() : main loop. Iterates the model forward and saves the state
            in the history files at every step
    """
    def __init__(self, param) :
        self.grid = grid.Grid(param) #loads the grid
        self.IO = nylesIO.NylesIO(param) #loads the IO
        self.initiate(param) #initiates the model, IO and needed variables


    def initiate(self,param) :
        if param['modelname'] == 'LES' :
            self.model = model_LES.LES(param)
        elif param['modelname'] == 'adv' :
            self.model = model_adv.Advection(param)

        self.tend = param['tend']
        self.cfl = param['cfl']
        self.dhist = param['timestep_history']

    def run(self) :
        print("Creating output file:", self.IO.hist_path)
        self.IO.init(self.model.state, self.grid)

        self.dx = 1. #Hardcoded for now TODO : get it from grid
        self.dt = self.cfl*self.dx

        tstep = np.arange(0,self.tend,self.dt)
        n = 0

        for t in tstep :
            n += 1
            self.model.forward(t, self.dt)
            #write to IO everytime t is multiple of dhist
            if (t%self.dhist == 0) :
                self.IO.do(self.model.state, t, n)
