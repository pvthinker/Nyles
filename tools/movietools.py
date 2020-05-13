import subprocess
import os

devnull = open(os.devnull, 'wb')

class Movie():
    """ Home made class to generate mp4 """

    def __init__(self, fig, name='mymovie', framerate=30):
        """ input: fig is the handle to the figure """
        self.fig = fig
        #dpi = fig.get_dpi()
        canvas_width, canvas_height = self.fig.canvas.get_width_height()
        #canvas_width *= dpi/100
        #canvas_height *= dpi/100
        print()
        print("dpi:", fig.get_dpi())
        print("fig:", fig)
        print("fig size inches:", fig.get_size_inches())
        print("canvas widthxheight:", fig.canvas.get_width_height())

        # Open an ffmpeg process
        outf = '%s.mp4' % name
        videoencoder = None
        for v in ['ffmpeg', 'avconv']:
            if subprocess.call(['which', v], stdout=subprocess.PIPE) == 0:
                videoencoder = v

        if videoencoder is None:
            print('\n')
            print('Neither avconv or ffmpeg was found')
            print('Install one or set param.generate_mp4 = False')
            raise ValueError('Install avconv or ffmpeg')
        else:
            print("Encoder is %s" % videoencoder)

        codec = "h264"#"mjpeg"#"mpeg4" #'libx264'
        
        check_codec_encoder(videoencoder, codec)
        
        cmdstring = (videoencoder,
                     '-y', '-r', str(framerate),  # overwrite, 30fps
                     # size of image string
                     '-s', '%dx%d' % (canvas_width, canvas_height),
                     '-pix_fmt', 'argb',  # format
                     '-f', 'rawvideo',
                     # tell ffmpeg to expect raw video from the pipe
                     '-i', '-',
                     '-vcodec', codec, outf)  # output encoding

        self.process = subprocess.Popen(cmdstring,
                                        stdin=subprocess.PIPE,
                                        stdout=devnull,
                                        stderr=devnull)

    def addframe(self):
        self.fig.set_dpi(100)
        print(self.fig.canvas.get_width_height())
        string = self.fig.canvas.tostring_argb()
        self.process.stdin.write(string)

    def finalize(self):
        self.process.communicate()

def check_codec_encoder(videoencoder, codec):
    res=subprocess.check_output([videoencoder, "-codecs"], stderr=subprocess.DEVNULL)
    lines=res.decode('utf-8').split('\n')
    ok = False
    for l in lines:
        s = l.split()
        if len(s) >=3:
            if codec in s[1]:
                if "E" == s[0][1]:
                    ok = True
    if not(ok):
        raise ValueError("%s cannot encode in %s" % (videoencoder, codec))
    else:
        print("encoding with codec '%s'" % codec)
    
if __name__ == '__main__':
    """

    Below is an example of how to use this module

    1/ prepare a figure (get the handles of it)
    2/ create the movie instance
    3/ do a loop and update all the handles
    4/ add each frame
    5/ finalize at then end of the loop

    """
    import numpy as np
    import matplotlib.pyplot as plt

    def f(t):
        """ the function to animate """
        return np.sin(2*np.pi*(x-t))

    # Create the plot
    x = np.linspace(-1, 1, 201)
    t = 0.
    template = 'time = %.2f'
    fig = plt.figure()
    fig.set_dpi(200)
    print(fig)
    print(fig.get_size_inches())
    print(fig.get_dpi())
    print(fig.canvas.get_width_height())
    p = plt.plot(x, f(t))
    plt.axis([-1, 1, -1, 1])
    p = p[0]
    plt.xlabel('x')
    plt.grid()
    ti = plt.title(template % t)

    # note that you don't even need to plot on the screen
    # to generate the movie!

    # create a Movie instance
    movie = Movie(fig, name='test_moving_sine')

    nt = 100
    dt = 1./nt
    for kt in range(nt):
        t = kt*dt
        # it's faster to update handles
        # than to recreate a plot!
        p.set_data(x, f(t))
        ti.set_text(template % t)
        fig.canvas.draw()
        # add the frame
        movie.addframe()

    # close the movie cleanly (mandatory)
    movie.finalize()
