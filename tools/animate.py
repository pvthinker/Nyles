#! /usr/bin/env python3

# The first line allows to run this program from the Linux shell without
# calling Python explicitly.  Since it is a comment, it is ignored on
# platforms which don't understand the "Hashbang".

"""Animate a Nyles history file and create a video from it.

This is a standalone Python3 script that can be called from the shell.
Call the script with the argument "--help" to see its usage.
"""

# Standard library imports
import os
import time
import argparse

# Third party imports
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib import animation

# Local imports
try:
    from plotting import Plotting
    from grid import Grid
    from variables import Scalar, State
except ModuleNotFoundError:
    print("Activate Nyles before running this script.")
    raise


## Settings for the video
#
# A common frame rate for artistic movies is 24 fps.  Use a lower frame
# rate if only few snapshots were saved.  Use a higher framerate
# (e.g. 30 or 60 fps) to have a more dynamic video if many snapshots
# were saved.
FPS = 4
# The standard aspect ratio for YouTube is 16:9 and a common resolution
# is 1080p (Full HD, 1920×1080 pixel).  For this resolution and a frame
# rate of 30 fps or less, YouTube recommends a bitrate of 8 Mbps.
# Reference: https://support.google.com/youtube/answer/1722171
RESOLUTION = [1920, 1080]
BPS = 8e6
# DPI determines how big the window of the animation is.  A bigger value
# will result in a smaller window and vice versa.  Thus increasing DPI
# is an easy way to have bigger fonts in the video and vice versa.
# Note: if the window is too small, it won't be possible to move the
# camera into a good angle; if the window is too big, it will be
# maximized and won't preserve the right aspect ratio.
DPI = 250
# Your name and institute to be saved in the metadata of the video.
CREATOR = ""


class Animator:
    """Creator for animations and videos."""

    def __init__(self, hist_path: str, varname="b", video_path=None, visible=True):
        """Open the history file and load parts of it."""
        # Initialize self.hist_file to prevent the destructor from failing
        self.hist_file = None

        # TODO: remove this when more variables are implemented
        if varname != "b":
            raise NotImplementedError("The plotting module does not yet support other variables than buoyancy.")

        # Save necessary arguments
        self.video_path = video_path
        self.visible = visible

        # Create the metadata for the video
        if self.video_path:
            # Extract the name of the experiment
            exp_name = os.path.basename(
                hist_path[:-8] if hist_path.endswith("_hist.nc") else hist_path
            )
            self.metadata = {
                "title": "Nyles experiment {}".format(exp_name),
                "artist": CREATOR,
                "genre": "Computational Fluid Dynamics (CFD)",
                "comment": "Created on {} with Nyles.  Nyles is a Large Eddy "
                           "Simulation written in Python.  For more information"
                           " visit https://github.com/pvthinker/Nyles."
                           .format(time.strftime('%d %b %Y')),
                "date": time.strftime("%Y-%m-%d"),
            }

        # Set parameters needed for Plotting
        # TODO: add more when plotting takes more parameters
        param = {
            "figsize": (RESOLUTION[0] / DPI, RESOLUTION[1] / DPI),
            "aspect": "equal",
        }

        # Open the history file and keep it open to allow sequential reading
        print("Loading history file {!r}:".format(hist_path))
        self.hist_file = nc.Dataset(hist_path)
        print(self.hist_file)

        # Load the needed data
        self.vardata = self.hist_file[varname]
        self.t = self.hist_file["t"]
        self.n = self.hist_file["n"]
        self.n_frames = self.n.size

        # Load parameters needed for Grid
        param["Lx"] = self.hist_file.Lx
        param["Ly"] = self.hist_file.Ly
        param["Lz"] = self.hist_file.Lz
        param["nx"] = self.hist_file.global_nx
        param["ny"] = self.hist_file.global_ny
        param["nz"] = self.hist_file.global_nz

        # Set parameters needed for Scalar
        param["nh"] = 0
        param["neighbours"] = {}

        # Create a Grid, the variable, a State, and a Plotting object
        grid = Grid(param)
        # Scalar takes actually a dimension instead of a unit, but this
        # does not matter as long as the units don't change
        scalar = Scalar(param, self.vardata.long_name, varname, self.vardata.units)
        state = State([scalar])
        self.p = Plotting(param, state, grid)

        # Get access to the 3D data array
        self.array = scalar.view("i")

    def __del__(self):
        """Close the history file in the destructor."""
        if self.hist_file:
            self.hist_file.close()

    def init(self):
        """Show the inital frame."""
        print("Variable:", self.vardata.long_name)
        print("Number of frames:", self.n_frames)
        if self.video_path:
            print("Output file:", self.video_path, end="")
            if os.path.exists(self.video_path):
                print(" -- file exists already and will be overwritten!")
            else:
                print("")
            if not self.visible:
                print("Fast mode: no animation will be visible during the process.")
            else:
                print('Slow mode: call script with "--fast" to speed up the video creation.')
        else:
            print("No video will be created.")
        # Load the initial data and show it
        self.array[...] = self.vardata[0]
        self.p.init(self.t[0], self.n[0])

    def run(self):
        """Create the animation and optionally save it."""
        if not self.video_path:
            plt.ioff()
        self.anim = animation.FuncAnimation(
            self.p.fig,
            self.update,
            frames=self.n_frames,
            repeat=False,
            interval=0,
        )
        if self.visible:
            plt.show()
        if self.video_path:
            self.anim.save(
                self.video_path,
                fps=FPS,
                dpi=DPI,
                bitrate=BPS,
                metadata=self.metadata,
            )

    def update(self, frame):
        """Load the data of the given frame and display it."""
        print("\rProcessing frame {} of {} ...".format(frame+1, self.n_frames), end="")
        # Load the data and show it
        self.array[...] = self.vardata[frame]
        self.p.update(self.t[frame], self.n[frame], self.visible)
        # At the end
        if frame + 1 == self.n_frames:
            if self.video_path:
                print("\b\b\b-- saved.")
            else:
                print("\b\b\b-- finished.")
                plt.pause(0.5)
                plt.close(self.p.fig)


if __name__ == "__main__":
    # Set up the command line arguments
    parser = argparse.ArgumentParser(
        description="Animate a Nyles history file and create a video from it."
    )
    parser.add_argument(
        dest="history_file", type=str,
        help="path to the Nyles history file",
    )
    parser.add_argument(
        "--out", dest="video_file", type=str,
        help="name of the video file to create; if this argument is omitted, "
             "the animation is shown, but no video is created",
    )
    parser.add_argument(
        "--fast", dest="fast", action="store_const", const=True, default=False,
        help="only save the video without showing it at the same time "
             "(this has no effect if no video file is given)",
    )
    parser.add_argument(
        "--var", dest="nickname", type=str,
        help='short name of the variable to animate (default: "b" for buoyancy)',
    )
    # TODO: make use of these arguments
    # parser.add_argument(
    #     "--cmap", dest="cmap", type=str,
    #     help="name of the colormap to use (default: coolwarm)",
    # )
    # parser.add_argument(
    #     "--vmax", dest="vmax", type=float,
    #     help="maximal value of the color axis (default: auto)",
    # )
    # parser.add_argument(
    #     "--vmin", dest="vmin", type=float,
    #     help="minimal value of the color axis (default: auto)",
    # )
    # Get the command line arguments
    args = parser.parse_args()
    hist_file_path = args.history_file
    video_file_path = args.video_file
    varname = args.nickname if args.nickname else "b"
    if video_file_path and args.fast:
        visible = False
    else:
        visible = True
    # TODO: make use of these arguments
    # cmap = args.cmap if args.cmap else "coolwarm"
    # vmax = args.vmax
    # vmin = args.vmin

    animator = Animator(hist_file_path, varname, video_file_path, visible)
    animator.init()

    print("")
    print("Do not close the figure!")
    print("Move the camera into a good angle,")
    input("then press Enter to start (or Ctrl+D to cancel). ")

    animator.run()
