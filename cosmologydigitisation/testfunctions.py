from numpy.fft import *
from matplotlib.pyplot import *
from numpy import *
from cosmdigitclasses import *
from matplotlib.image import *

A = normal(size = 100)
powa = sum(A*A)
print("A: " + str(powa))

B = digitise1bit(A, 5)
powb = sum(B*B)
print("B: " + str(powb))
B = (powa/powb)**0.5 * B
powb = sum(B*B)
print("B: " + str(powb))