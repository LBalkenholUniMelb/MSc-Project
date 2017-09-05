from numpy.fft import *
from matplotlib.pyplot import *
from numpy import *
from cosmdigitclasses import *
from matplotlib.image import *

A = normal(size = 100)

plot(range(100), A)
plot(range(100), digitise1bit(A, 1))
show()