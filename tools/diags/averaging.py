import numpy as np

def horizontal_avg(field, yidx=slice(None), xidx=slice(None)):
    """
    Do a horizontal average

    Parameters
    ----------
    field : array,  2D or 3D
         the field that you want to average

    xidx, yidx : slice, optional
         ranging over the interior points, excluding the halo
         default is all points


    Returns
    -------
    res : float or 1D array
       the horizontal average

    """
    if len(field.shape) == 2:
        res = np.mean(field[yidx, xidx])
    elif len(field.shape) == 3:
        res = np.mean(np.mean(field[:, yidx, xidx], axis=2), axis=1)
    else:
        raise ValueError("field should be a 2D or 3D array")
    return res


def azimuthal_avg(field, x, y, dx, dy):

    # todo: make the function accepts 3D or 4D array for field
    assert len(field.shape) == 2
    assert len(x.shape) == 2
    assert len(y.shape) == 2
    
    d = np.sqrt(x**2+y**2)

    
    xmax = np.max(x)
    ymax = np.max(y)
    rmax = min([xmax, ymax])

    dr = min(dx, dy)
    nr = int(rmax/dr+0.5)
    
    r = np.arange(nr)*dr
    # mid ring radial distance
    rm = 0.5*(r[1:]+r[:-1])
    
    res = np.zeros(nr-1)
    # average over rings
    for k in range(nr-1):
        res[k] = np.mean(field[(d >= r[k]) & (d < r[k+1])])

    return rm, res

def test_azimuthal_avg():
    nx, ny = 201, 201
    x = np.linspace(-1,1,nx)
    y = np.linspace(-1,1,ny)

    xx, yy = np.meshgrid(x, y)

    r = np.sqrt(xx**2+yy**2)

    f = lambda x: np.sin(x*5)
    
    phi = f(r)
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    rm, phiavg = azimuthal_avg(phi, xx, yy, dx, dy)

    plt.clf()
    plt.plot(rm, phiavg, '+', label='computed')
    plt.plot(rm, f(rm), label='exact')
    plt.legend()
    plt.xlabel('r')
    
if __name__ == '__main__':

    import matplotlib.pyplot as plt

    plt.ion()
    
    test_azimuthal_avg()
