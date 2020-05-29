from __future__ import print_function
import matplotlib as mpl

font = {'size': 14}
mpl.rc('font', **font)
mpl.rcParams['toolbar'] = 'None'
mpl.use('TkAgg') # <= the backend is crucial
# figure.canvas.to_string_rgb() should return 24 bits image
# the default Qt5Agg doesn't work

from shutil import which
import argparse
import glob
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
from iosubdomains import Variable
import movietools

def init_figure(fig, str_time, field2d, Lx, Ly):
    """create the figure

    this is where you adapt the figure to your needs
    the function should return the graphical objects to update
    """

    field_size = np.shape(field2d)
    zoom_x = (0.8*fig_size[0])//field_size[1]
    zoom_y = (0.9*fig_size[1])//field_size[0]
    # print("zoom factors x : %i / y : %i" % (zoom_x, zoom_y))
    zoom_factor = min(4, min(zoom_x, zoom_y))
    print("each grid cell is (%i, %i) pixels" % (zoom_factor, zoom_factor))
    rectangle = [0.1, 0.1,
                 field_size[1]/fig_size[0]*zoom_factor,
                 field_size[0]/fig_size[1]*zoom_factor]
    ax = plt.axes()#rectangle)

    im = ax.imshow(field2d, cmap=plt.get_cmap('RdBu_r'),
                   vmin=cax[0], vmax=cax[1], origin='lower',
                   extent=[0, Lx, 0, Ly])
    ti = ax.set_title(str_time)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    #plt.tight_layout()
    #divider = make_axes_locatable(ax)
    #cbax = divider.append_axes("right", size="3%", pad=0.1)

#    pos = [rectangle[0]+rectangle[2]+0.02, rectangle[1], 0.05, rectangle[3]]
    cb = fig.colorbar(im)

    #pos = np.array(cb.ax.get_position().bounds)
    #pos[1], pos[3] = rectangle[1], rectangle[3]
    #plt.tight_layout()
#    cb.ax.set_position(pos)
#    fig.canvas.draw()
    return im, ti

def time_into_string(time):
    """ convert the time into a suitable string for the title"""

    t = time/86400.
    day = int(t)
    t1 = (t-day)*24
    hour = int(t1)
    t2 = (t1-hour)*60
    minute = int(t2)
    string = '%id-%i:%02i' % (day, hour, minute)
    string = '%05.3f' % time
    return string

if __name__ == '__main__':

    # https://docs.python.org/2/library/argparse.html#module-argparse
    parser = argparse.ArgumentParser(description='Do a movie from netfiles')
    parser.add_argument('--output', dest='video_file', type=str,
                        help='name of the video file')
    parser.add_argument('--var', dest='varname', type=str, required=True,
                        help='variable to animate')
    parser.add_argument('--cax', dest='cax', type=float, nargs=2,
                        help='two floats for the color range')
    parser.add_argument('--kz', dest='kz', type=int, required=True,
                        help='level to plot')
    parser.add_argument('--dpi', dest='dpi', type=int,
                        help='dpi of the image (default 100)')
    parser.add_argument('--dir', dest='direxp', type=str, required=True,
                        help='directory of the netCDF files')
    parser.add_argument('-auto', dest='auto', action="store_true",
                        help='autoscale for color range')

    # save the user options
    args = parser.parse_args()
    varname = args.varname
    cax = args.cax
    direxp = os.path.join(args.direxp, '')
    if cax is None:
        cax = (None, None)
    kz = args.kz
    if args.dpi is None:
        my_dpi = 100
    else:
        my_dpi = args.dpi
    if args.video_file is None:
        video_file = '%s.mp4' % varname
    else:
        video_file = args.video_file
    autoscale = args.auto
    # if args.auto is None:
    #     autoscale = False
    # else:
    #     autoscale = True
    # template string for the figure title
    title_string = '%s / time = '

    # control the options
    print()
    print('varname = %s' % varname)
    print('level = %i' % kz)
    print('cax = ', cax)
    print("directory = ", direxp)
    print()

    plt.ion()
    zoom_factor = 1
    # best youtube aspect ratio is 16:9
    fig_size = np.array([1280, 720])
    #figsize = (12.8, 7.2)
    figsize = (6.4, 4.8)
    fig = plt.figure()
    fig.set_dpi(my_dpi)
    fig.set_size_inches(figsize)

    nc00 = glob.glob(direxp+"*00_hist.nc")[0]
    ncall = glob.glob(direxp+"*_hist.nc")

    nctemplate = nc00.replace("00_hist", "%02i_hist")

    #fig = plt.figure(figsize=fig_size/my_dpi, dpi=my_dpi)

    time = Variable(nctemplate, "t")[:]
    with Dataset(nc00) as nc:
        Lx = nc.getncattr("Lx")
        Ly = nc.getncattr("Ly")
    V = Variable(nctemplate, varname)
    nb_frames = V.shape[0]
    print("Template    : %s" % nctemplate)
    print("Nb of frames: %i" % nb_frames)
    print("Nb of files : %i" % len(ncall))

    if autoscale:
        mini = 9e99
        maxi = -9e99
        for kt in range(0, nb_frames, int(nb_frames/4)):
            field2d = V[kt,kz]
            mini = min(mini, np.min(field2d))
            maxi = max(maxi, np.max(field2d))
        delta = maxi-mini
        maxi += delta/5
        mini -= delta/5
        print("Min - Max : %.2g - %.2g" % (mini, maxi))
        cax = [mini, maxi]
    #plt.show(block=False)    
    #plt.savefig("z.png")

    depth = "z = %.0f m" % Variable(nctemplate, "z")[kz]
    template = " / ".join([varname, depth, "time = %4.2f day"])
    day = 86400.
    
    for kt in range(nb_frames):
        field2d = V[kt,kz]
        t_str = template % (time[kt]/day)
        if kt == 0:
            im, ti = init_figure(fig, t_str, field2d, Lx, Ly)
            movie = movietools.Movie(fig, "movie")
        else:
            im.set_data(field2d)
            ti.set_text(t_str)
        fig.canvas.draw()

        movie.addframe()
    print()

    movie.finalize()

        
