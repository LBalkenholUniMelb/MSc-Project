from numpy.fft import *
from matplotlib.pyplot import *
from numpy import *
from cosmdigitclasses import *
from matplotlib.image import *

pn = 1024
cmbnoise = normal(size = (pn, pn*2))
#savetxt("cmbnoiseseed.txt", cmbnoise)

# Read in data
file = open("cmbnoiseseed.txt")
re = zeros((pn, pn))
im = zeros((pn, pn))
rowi = 0
for row in file:
    rowvalues = [float(i) for i in row.split()]
    re[rowi] = rowvalues[:pn]
    im[rowi] = rowvalues[pn:]
    rowi += 1

print(shape(re))
print(re)
print(shape(im))
print(im)