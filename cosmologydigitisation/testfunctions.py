from numpy.fft import *
from matplotlib.pyplot import *
from numpy import *
from numpy.random import *
from copy import *
from cosmdigitclasses import *
from matplotlib.image import *
from digitisationschemes import *
from csv import *

pixelnumber = 512
file = open("CMBMap.txt")
cmbmap = zeros((pixelnumber, pixelnumber))
rowindex = 0
for row in file:
    rowvalues = [float(i) for i in row.split()]
    cmbmap[rowindex] = rowvalues[:pixelnumber]
    rowindex += 1



wr = writer(open("CMBMapascsv.csv", 'w'))
for row in cmbmap:
    wr.writerow(row)