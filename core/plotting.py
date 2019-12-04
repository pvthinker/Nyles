from enum import Enum

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from timing import timing


class PlotType(Enum):
    B_Interface = 0
    Tracer = 1


class Plotting(object):

    COLORS = ["b", "r", "g", "y", "c", "m", "k"]

    def __init__(self, param, state, grid):
        self.state = state
        self.grid = grid

        # Parameter figsize is not used for live-plotting, only for the
        # creation of videos in post-processing
        self.figsize = param["figsize"] if "figsize" in param else None
        # Parameter n_update is not used in the animate module because
        # the frames saved in the history file might not be multiples of
        # any integer
        self.n_update = (
            param["iterations_per_frame"] if "iterations_per_frame" in param else 1
        )
        self.aspect = param["aspect"]
        self.rotation_speed = param["rotation_speed"]

        if param["style"] == "b-interface":
            self.plot_type = PlotType.B_Interface
            self.variable = self.state.b
            self.data = self.variable.view("i")
            self.show_bigger = param["stable_stratification"]
            # For a surface, 2D coordinates are needed
            self.x_2D = self.grid.x_b.view("i")[0, :, :]
            self.y_2D = self.grid.y_b.view("i")[0, :, :]
            self.z_2D = np.zeros_like(self.x_2D)
            self.z_3D = self.grid.z_b.view("i")[:, :, :]
        elif param["style"] == "tracer":
            self.plot_type = PlotType.Tracer
            self.n_tracers = param["n_tracers"]
            self.variables = [self.state.get("t{}".format(i)) for i in range(self.n_tracers)]
            # For a scatter plot, 3D coordinates are needed
            self.x_3D = self.grid.x_b.view("i")[:, :, :]
            self.y_3D = self.grid.y_b.view("i")[:, :, :]
            self.z_3D = self.grid.z_b.view("i")[:, :, :]
            # Pre-allocate arrays
            self.condition = np.zeros_like(self.x_3D)
            self.x_coords = [np.zeros_like(self.x_3D) for i in range(self.n_tracers)]
            self.y_coords = [np.zeros_like(self.y_3D) for i in range(self.n_tracers)]
            self.z_coords = [np.zeros_like(self.z_3D) for i in range(self.n_tracers)]

    def init(self, t, n):
        if self.plot_type == PlotType.B_Interface:
            # Get central value of the data; this is where the interface is shown
            self.data[...] = self.variable.view("i")
            self.level = (np.min(self.data) + np.max(self.data)) / 2

        plt.ion()
        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.gca(projection="3d", azim=-135, elev=10)
        self.ax.set_title("n = {}, t = {:.2f}".format(n, t), fontdict={"family": "monospace"})
        self.update(t, n)

    def update(self, t, n, visible=True):
        if n % self.n_update:
            # Only update if n is a multiple of self.n_update
            return
        if self.plot_type == PlotType.B_Interface:
            self.calculate_interface()
        elif self.plot_type == PlotType.Tracer:
            self.calculate_coordinates()
        self.plot_data(t, n, visible)

    @timing
    def plot_data(self, t, n, visible=True):
        self.ax.clear()

        self.ax.title.set_text("n = {}, t = {:.2f}".format(n, t))
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_zlabel("z")

        self.ax.set_aspect(self.aspect)
        self.ax.set_xlim(0, self.grid.Lx)
        self.ax.set_ylim(0, self.grid.Ly)
        self.ax.set_zlim(0, self.grid.Lz)

        # Rotate the camera around the 3D scene
        self.ax.azim += self.rotation_speed
        # Other camera settings can be configured with self.ax.elev and
        # self.ax.dist, for example:
        # self.p.ax.elev = 10 * np.sin(2*np.pi*frame/60)
        # self.p.ax.dist = 10 + 5 * np.sin(2*np.pi*frame/120)

        if self.plot_type == PlotType.B_Interface:
            self.ax.plot_surface(self.x_2D, self.y_2D, self.z_2D)
            # Draw contour diagrams on the correct faces
            self.ax.azim %= 360
            if self.ax.azim >= 270:
                xoff = 0
                yoff = self.grid.Ly
            elif self.ax.azim > 180:
                xoff = self.grid.Lx
                yoff = self.grid.Ly
            elif self.ax.azim > 90:
                xoff = self.grid.Lx
                yoff = 0
            elif self.ax.azim >= 0:
                xoff = 0
                yoff = 0
            self.ax.contourf(
                self.x_2D, self.y_2D, self.z_2D, zdir='z', offset=0, cmap="coolwarm"
            )
            self.ax.contourf(
                self.x_2D, self.y_2D, self.z_2D, zdir='x', offset=xoff, cmap="coolwarm"
            )
            self.ax.contourf(
                self.x_2D, self.y_2D, self.z_2D, zdir='y', offset=yoff, cmap="coolwarm"
            )
        elif self.plot_type == PlotType.Tracer:
            for t in range(self.n_tracers):
                self.ax.scatter(
                    self.x_coords[t], self.y_coords[t], self.z_coords[t],
                    c=self.COLORS[t%7],
                )

        if visible:
            self.fig.canvas.draw()

    @timing
    def calculate_coordinates(self):
        for t in range(self.n_tracers):
            n = 1  # increase to make plot less dense
            condition = self.variables[t].view("i")[::n, ::n, ::n] > 0.5
            self.x_coords[t] = self.x_3D[::n, ::n, ::n][condition]
            self.y_coords[t] = self.y_3D[::n, ::n, ::n][condition]
            self.z_coords[t] = self.z_3D[::n, ::n, ::n][condition]

    @timing
    def calculate_interface(self):
        self.data[...] = self.variable.view("i")
        for i in range(self.z_2D.shape[1]):
            for j in range(self.z_2D.shape[0]):
                indices = np.where(self.data[:, j, i] >= self.level)
                assert len(indices) == 1
                indices = indices[0]
                if indices.size > 0:
                    if self.show_bigger:
                        self.z_2D[j, i] = self.z_3D[indices[0], 0, 0]
                    else:
                        self.z_2D[j, i] = self.z_3D[indices[-1], 0, 0]
                else:
                    self.z_2D[j, i] = self.z_3D[-1, 0, 0]
