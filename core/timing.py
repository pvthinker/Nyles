from functools import wraps
from time import time
import pickle
import numpy as np
import mpitools
import socket

# https://stackoverflow.com/questions/1622943/timeit-versus-timing-decorator
stats = {}

# list of hostnames where doing a plot fails
blacklist = ["irene"]


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        if "forward" in stats.keys():
            if len(stats["forward"]) >= 100:
                ok = True
            else:
                ok = False
        else:
            ok = False
        if ok:
            result = f(*args, **kw)
        else:
            ts = time()
            result = f(*args, **kw)
            te = time()
            if 'timing' in kw.keys():
                pass
            else:
                name = f.__name__
                if name in stats.keys():
                    stats[name] += [te-ts]
                else:
                    stats[name] = [te-ts]

            #print('func:%r took: %2.4e sec' % (f.__name__,  te-ts))
        return result
    return wrap

def write_timings(path):
    if mpitools.get_myrank() == 0:
        fid = open('%s/timing.pkl' % path, 'bw')
        pickle.dump(stats, fid)
    
def analyze_timing(path):
    host = socket.gethostname()
    skip = any([machine in host for machine in blacklist])

    if (mpitools.get_myrank() == 0) and not(skip):

        import matplotlib as mpl
        import matplotlib.pyplot as plt

        mpl.rcParams['font.size'] = 14
        mpl.rcParams['lines.linewidth'] = 2

        filename = '%s/timing.pkl' % path
        pngtiming = '%s/timing.png' % path

        f = open(filename, 'br')
        timing = pickle.load(f)

        mean = []
        keys = []
        for k, vals in timing.items():
            mean += [np.mean(vals)]
            keys += [k]

        idx = np.argsort(mean)

        plt.figure(figsize=(10, 5))
        for k in idx[::-1]:
            vals = timing[keys[k]]
            plt.loglog(vals, label=keys[k], alpha=0.8)

        plt.xlabel('iteration')
        plt.ylabel('time [s]')
        gca = plt.gca()
        gca.set_ylim([1e-5, 10])
        plt.grid()
        plt.legend(loc='upper right', fontsize=10)
        plt.savefig(pngtiming)
