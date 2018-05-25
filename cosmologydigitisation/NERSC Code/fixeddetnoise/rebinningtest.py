#--- Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from cosmdigitclasses import *
from numpy import *
from scipy.signal import convolve2d
from scipy.stats import *
from pickle import *
from scipy.optimize import curve_fit, leastsq

rc('text', usetex=True)
rc("xtick", labelsize = 20)
rc("ytick", labelsize = 20)

cosm = Cosmologist()


realmap = normal(0, 1, size = (512, 512))


k, p, err, h = cosm.sriniPowerSpectrum([512, 512, 2, 2], realmap)
#ke, pe, erre, he = cosm.sriniPowerSpectrum([512, 512, 2, 2], realmap)

plot(k, p)
xscale("linear")
yscale("linear")
show()