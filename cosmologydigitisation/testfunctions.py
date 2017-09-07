from numpy.fft import *
from matplotlib.pyplot import *
from numpy import *
from cosmdigitclasses import *
from matplotlib.image import *


A = normal(size = 10)

A = digitise2bithalfmax(A, 0)

plot(range(10), A)
show()