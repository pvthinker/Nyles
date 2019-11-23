from functools import wraps
from time import time

# https://stackoverflow.com/questions/1622943/timeit-versus-timing-decorator

def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        #print('func:%r took: %2.4e sec' % (f.__name__,  te-ts))
        return result
    return wrap

