import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from timing import timing


class Plotting(object):
    def __init__(self, param, state, grid):
        # Parameters for the video creation
        self.figsize = param["figsize"] if "figsize" in param else None
        self.aspect = param["aspect"]
        # TODO: make more use of param

        self.state = state
        self.grid = grid

        self.plot_variable = self.state.b

        # For variables at b-point (not velocity, not vorticity)
        self.X = self.grid.x_b.view("i")[0, :, :]
        self.Y = self.grid.y_b.view("i")[0, :, :]
        self.Z = np.zeros_like(self.X)
        self.z_3D = self.grid.z_b.view("i")[:, :, :]

        # Put the following parameter to True for situations with stable
        # stratification, use False for unstable stratifications
        self.show_bigger = False

    def init(self, t, n):
        self.data_array = self.plot_variable.view("i")
        self.level = (np.min(self.data_array) + np.max(self.data_array)) / 2

        plt.ion()
        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.add_subplot(1, 1, 1, projection="3d")
        self.ax.set_title("n = {}, t = {:.2f}".format(n, t), fontdict={"family": "monospace"})
        self.update(t, n)

    def update(self, t, n, visible=True):
        self.calc_data()
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

        self.ax.plot_surface(self.X, self.Y, self.Z)
        self.ax.contourf(self.X, self.Y, self.Z, zdir='z', offset=0,
                         cmap="coolwarm")
        self.ax.contourf(self.X, self.Y, self.Z, zdir='x', offset=0,
                         cmap="coolwarm")
        self.ax.contourf(self.X, self.Y, self.Z, zdir='y', offset=self.grid.Ly,
                         cmap="coolwarm")

        if visible:
            self.fig.canvas.draw()

    @timing
    def calc_data(self):
        self.data_array = self.plot_variable.view("i")
        for i in range(self.Z.shape[1]):
            for j in range(self.Z.shape[0]):
                indices = np.where(self.data_array[:, j, i] >= self.level)
                assert len(indices) == 1
                indices = indices[0]
                if indices.size > 0:
                    if self.show_bigger:
                        self.Z[j, i] = self.z_3D[indices[0], 0, 0]
                    else:
                        self.Z[j, i] = self.z_3D[indices[-1], 0, 0]
                else:
                    self.Z[j, i] = self.z_3D[-1, 0, 0]
