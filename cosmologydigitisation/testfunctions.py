from numpy.fft import *
from matplotlib.pyplot import *
from numpy import *
from numpy.random import *
from copy import *
from cosmdigitclasses import *
from matplotlib.image import *
from digitisationschemes import *

Ar = normal(size = (100, 100))
Ai = normal(size = (100, 100))
freqx = fftshift(fftfreq(100, 1))
freqx2d = tile(freqx, (100, 100))
freqy = transpose(freqx)
freqy2d = transpose(tile(freqy, (100, 100)))

cl2d = zeros((100, 100))
lknee = 2
for y in range(100):
    for x in range(100):
        l = (2*pi*freqx2d[y][x]**2 + 2*pi*freqy2d[y][x]**2)**0.5
        cl2d[y][x] = (1 + (l/lknee)**(8.0/3.0))
    print(y)
print(freqx)
noisenorm = real(ifft2(Ar + 1j*Ai))

Ar = Ar * cl2d
Ai = Ai * cl2d

noise = real(ifft2(Ar + 1j*Ai))

c = Cosmologist()
k, p, perr, hits = c.sriniPowerSpectrum([100, 100, 1, 1], noise)
k0, p0, perr0, hits0 = c.sriniPowerSpectrum([100, 100, 1, 1], noisenorm)
plot(k, p, "r")
plot(k0, p0, "b")
show()