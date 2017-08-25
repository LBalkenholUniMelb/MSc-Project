from numpy.fft import *
from matplotlib.pyplot import *
from numpy import *
from cosmdigitclasses import *

A = ones((64, 64))
Ap = zeropad(A)
imshow(Ap)
colorbar()
show()