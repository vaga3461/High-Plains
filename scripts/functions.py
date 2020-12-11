""" This function does down-to-the-east tilting """
def down_to_east(x_z, scalar):    
    return (x_z*scalar)

""" This functions does down-to-the-northest tilting """
def down_to_northeast(x, y, x_scale, y_scale):
    xresult = x*x_scale
    yresult = y*y_scale
    return xresult+yresult

""" This functions does gaussian tilting
Must import numpy and scipy.stats as stats to use """
def gaussian(x_raster, variance, sig_scale, maximum):
    mu = x_raster[:]
    sigma = np.sqrt(variance)
    x = np.linspace(mu + sig_scale*sigma, mu + 0*sigma, maximum)
    curve = maximum*stats.norm.pdf(x, mu, sigma)
    return curve

""" This is a test """
def test_add(x, y):
    return x+y
    