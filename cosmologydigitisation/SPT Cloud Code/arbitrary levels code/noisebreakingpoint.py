# #--- Setup MPI
# from mpi4py import MPI
# comm = MPI.COMM_WORLD
# processorrank = comm.Get_rank()

#--- Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from scipy.stats import binned_statistic
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from cosmdigitclasses import *
from numpy import *
from scipy.signal import convolve2d
from scipy.stats import *


cosm = Cosmologist()

#--- Define map parameters
fieldsizearcmins = 2048
pixelsizearcmin = 2
pixelnumber = int(fieldsizearcmins/pixelsizearcmin)
df = 1.0
mapfieldsize = int(fieldsizearcmins/2.0)
mappixelnumber = int(pixelnumber/2.0)
declims = [0, mapfieldsize] #arcmins
ralims = [0, mapfieldsize] #arcmins
readoutfreq = 200 #Hz
raspeed = 0.2 #arcmin/s
nodecscans = mappixelnumber
norablocks = mappixelnumber
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
observationlim = 1
pixelspacingrad = float(pixelsizearcmin)*arcmins2radians

#--- Define noise level
noisemap = 0.5 #3.0 # muK arcmin
noisepix = noisemap/float(pixelsizearcmin)
fsky = (mapfieldsize*mapfieldsize)/(4.0*pi*60.0*60.0*(180.0/pi)**2.0)
noisecl = 4.0*pi*fsky*noisepix*noisepix/float(mappixelnumber*mappixelnumber)
eta = sqrt(float(compression)) * sqrt(float(observationlim)) * sqrt(float(noisecl) / (float(pixelnumber * pixelnumber)))


print(compression)


#--- Read in CMP
cmbmap = zeros((mappixelnumber, mappixelnumber))
y = 0
for row in open("cmbmaparblvltemplate.txt"):
    rowvalues = row.split()
    for x in range(mappixelnumber):
        cmbmap[y, x] = float(rowvalues[x])
    y += 1

timestream = cmbmap.flatten()
timestream = repeat(timestream, compression)
timestream = timestream[0:norablocks*compression]
timestreamcl = fft.fft(timestream)
timestreamcl = timestreamcl*conjugate(timestreamcl)/(float(len(timestream)))
timestreaml = 2.0*pi*fft.fftfreq(len(timestream), float(pixelspacingrad)/float(compression))

#timestreamnoise = normal(0, eta, len(timestream))
#timestreamnoisecl = fft.fft(timestreamnoise)
#timestreamnoisecl = timestreamnoisecl*conjugate(timestreamnoisecl)/(float(len(timestream)))


plot(timestreaml, timestreamcl, "b")
#plot(timestreaml, timestreamnoisecl, "r")
axvline(x = 5400, color = "k")
axhline(y = eta**2.0, color = "r")
title("CMB timestream PSD Simplified")
#xlim((20, 6500))
#ylim((1e-17, 1e-5))
xscale("log")
yscale("log")
show()