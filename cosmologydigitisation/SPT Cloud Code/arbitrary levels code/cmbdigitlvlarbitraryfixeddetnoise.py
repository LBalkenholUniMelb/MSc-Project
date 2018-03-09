# #--- Setup MPI
#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#processorrank = comm.Get_rank()

#--- Define digitisation functions

def digitise1bit(signal, lvl):
    comp = signal < 0
    signal[:] = lvl
    signal[comp] = -1.0*lvl
    return signal

#--- Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from scipy.stats import binned_statistic
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from numpy.random import *
from numpy import *
from scipy.signal import convolve2d
from scipy.stats import *


arcmins2radians = np.radians(1./60.)

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
raspeed = 0.0005 #0.0005 #arcmin/s
nodecscans = mappixelnumber
norablocks = mappixelnumber
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
observationlim = 1


cmbmap = zeros((mappixelnumber, mappixelnumber))
fname = "cmbtemplatemapfixeddetnoise.txt"
y = 0
for row in open(fname):
    rowvalues = row.split()
    for x in range(mappixelnumber):
        cmbmap[y, x] = float(rowvalues[x])
    y += 1


#--- Define Noise

noisedet = 500.0 # muK sqrt(s)
noiseinduced = noisedet/(sqrt(1.0/float(readoutfreq))) # muK
pixelrms = noisedet/sqrt(float(pixelsizearcmin)/raspeed) # muK
noisemap = pixelrms*float(pixelsizearcmin) # muK arcmin
noisecl = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixelrms*pixelrms # muK^2


#--- Recreate Observation
cmbnoisemap = zeros((nodecscans, norablocks))
cmbnoisemap1bit = zeros((nodecscans, norablocks))
#cmbnoisemap2bithm = zeros((nodecscans, norablocks))
#cmbnoisemap2bitopt = zeros((nodecscans, norablocks))


observationind = 0

while observationind < observationlim:

    # Decompose into bolometer signals
    for d in range(nodecscans):

        for ri in range(norablocks):

            # create noise
            noiser = normal(0, noiseinduced, compression)

            # add noise
            tod = cmbmap[d, ri] + noiser
            cmbnoisemap[d, ri] += mean(tod)

            # digitise
            tod1bit = digitise1bit(tod, 1.0)
            cmbnoisemap1bit[d, ri] += mean(tod1bit)


            #tod2bitopt = deepcopy(tod)

            #tod2bithm = digitise2bithalfmax(tod)
            #cmbnoisemap2bithm[d, ri] += mean(tod2bithm)

            #tod2bitopt = digitise2bitoptimal(tod2bitopt, eta)
            #cmbnoisemap2bitopt[d, ri] += mean(tod2bitopt)


        print("Noise, Digitisation and Compression: " + str(int(100*d/nodecscans)) + "%")

    print("--- COMPLETED OBSERVATION " + str(observationind) + " ---")

    observationind += 1

#--- Normalise Maps

cmbnoisemap = cmbnoisemap * (1.0/float(observationlim))
cmbnoisemap1bit = cmbnoisemap1bit * (1.0/float(observationlim))
#cmbnoisemap2bithm = cmbnoisemap2bithm * (1.0/float(observationlim))
#cmbnoisemap2bitopt = cmbnoisemap2bitopt * (1.0/float(observationlim))

#--- Save results

#savetxt("cmbmap" + str(processorrank) + ".txt", cmbmap)
savetxt("cmbnoisemapfixeddetnoise" + str(processorrank) + ".txt", cmbnoisemap)
savetxt("cmbnoisemap1bitfixeddetnoise" + str(processorrank) + ".txt", cmbnoisemap1bit)
#savetxt("cmbnoisemap2bithm" + str(processorrank+10) + ".txt", cmbnoisemap2bithm)
#savetxt("cmbnoisemap2bitopt" + str(processorrank+10) + ".txt", cmbnoisemap2bitopt)