from numpy.fft import *
from matplotlib.pyplot import *
from numpy import *
from cosmdigitclasses import *
from matplotlib.image import *

A = ones((64, 64))
B = zeros((64, 64))

fig = figure()
a = fig.add_subplot(1,2,1)
imgplot = imshow(A)
a.set_title('Before')
colorbar(orientation ='horizontal')
a = fig.add_subplot(1,2,2)
imgplot = imshow(B)
imgplot.set_clim(0.0,0.7)
a.set_title('After')
colorbar(orientation='horizontal')
show()