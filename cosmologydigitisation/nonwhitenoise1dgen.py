from numpy import *
from numpy.random import *
from numpy.fft import *
from scipy import *
from matplotlib.pyplot import *

def nonwhitenoise(l):
    lknee = 0.05
    if l == 0:
        return 1.0
    else:
        return (1.0 + (abs(l)/lknee)**(-8.0/3.0))#-8/3

scanhits = 100
hitspacing = 1

l = 2.0*pi*fftshift(fftfreq(scanhits, hitspacing))

print(max(l))

msk = asarray([nonwhitenoise(k) for k in l])

#plot(l, msk)
#show()

noisef = normal(0, 1, scanhits)

noisef = noisef*msk

noiser = real(ifft(noisef))

plot(range(scanhits), noiser)
show()