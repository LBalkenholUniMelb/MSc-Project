#--- SETUP TIME PROFILING
from time import *

#--- TIME
tsetupstart = time()

#--- Setup MPI
#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#processorrank = comm.Get_rank()

#--- Make necessary imports
from numpy import zeros, arange, real, shape, round
from numpy import *
from numpy.random import *
from copy import deepcopy
from pickle import *

from matplotlib.pyplot import *

#--- Define digitisation functions and helpers
def digitise1bitefficient(signal, lvl):
    return mean(((signal > 0)-0.5)*2*lvl)

def digitise2bitoptimalefficient(signal, var):
    return (1.0/float(len(signal))) * ( sum(signal >= 0.9816*var)*1.51*var + sum((signal >= 0) & (signal < 0.9816*var))*0.4528*var - sum(signal < -0.9816*var)*1.51*var - sum((signal < 0) & (signal >= -0.9816*var))*0.4528*var)

def digitise3bitoptimalefficient(signal, var):
    pos = (1.0/float(len(signal))) * ( sum(signal >= 1.748*var)*2.152*var + sum((signal >= 1.050*var) & (signal < 1.748*var))*1.344*var + sum((signal >= 0.501*var) & (signal < 1.050*var))*0.756*var + sum((signal >= 0) & (signal < 0.501*var))*0.245*var)
    neg = (1.0/float(len(signal))) * ( sum(signal < -1.748*var)*2.152*var + sum((signal < -1.050*var) & (signal >= -1.748*var))*1.344*var + sum((signal < -0.501*var) & (signal >= -1.050*var))*0.756*var + sum((signal < 0) & (signal >= -0.501*var))*0.245*var)
    return pos-neg

arcmins2radians = radians(1./60.)

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
raspeed = 0.0005 #arcmin/s
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


#--- TIME
tsetupstop = time()
tsetup = tsetupstop-tsetupstart
#print("----------" + str(processorrank) + "----------")
print("SETUP COMPLETED IN:")
print(tsetup)
#print("----------" + str(processorrank) + "----------")
tobsstart = time()


#--- Recreate Observation
cmbnoisemap = zeros((nodecscans, norablocks))
cmbnoisemap1bit = zeros((nodecscans, norablocks))
cmbnoisemap2bitopt = zeros((nodecscans, norablocks))
cmbnoisemap3bitopt = zeros((nodecscans, norablocks))


observationind = 0

while observationind < observationlim:

    # Decompose into bolometer signals
    for d in range(1):

        for ri in range(norablocks):

            # create noise
            noiser = normal(0, noiseinduced, compression)

            # add noise
            tod = cmbmap[d, ri] + noiser
            cmbnoisemap[d, ri] += mean(tod)

            # efficiently digitise
            cmbnoisemap1bit[d, ri] +=  digitise1bitefficient(tod, 1.0)
            cmbnoisemap2bitopt[d, ri] += digitise2bitoptimalefficient(tod, noiseinduced)
            cmbnoisemap3bitopt[d, ri] += digitise3bitoptimalefficient(tod, noiseinduced)

        print("Noise, Digitisation and Compression: " + str(int(100*d/nodecscans)) + "%")

    print("--- COMPLETED OBSERVATION " + str(observationind) + " ---")

    observationind += 1

#--- TIME
tobsstop = time()
tobs = tobsstop-tobsstart
#print("----------" + str(processorrank) + "----------")
print("1 ROW OBSERVATION COMPLETED IN:")
print(tobs)
#print("----------" + str(processorrank) + "----------")
tsavestart = time()

#--- Normalise Maps

cmbnoisemap = cmbnoisemap * (1.0/float(observationlim))
cmbnoisemap1bit = cmbnoisemap1bit * (1.0/float(observationlim))
cmbnoisemap2bitopt = cmbnoisemap2bitopt * (1.0/float(observationlim))
cmbnoisemap3bitopt = cmbnoisemap3bitopt * (1.0/float(observationlim))

#--- Save results

#dump(cmbnoisemap, open("cmbnoisemap" + str(processorrank) + ".p", "wb"), -1)
#dump(cmbnoisemap, open("cmbnoisemap1bit" + str(processorrank) + ".p", "wb"), -1)
#dump(cmbnoisemap, open("cmbnoisemap2bitopt" + str(processorrank) + ".p", "wb"), -1)
#dump(cmbnoisemap, open("cmbnoisemap3bitopt" + str(processorrank) + ".p", "wb"), -1)

#--- TIME
tsavestop = time()
tsave = tsavestop-tsavestart
#print("----------" + str(processorrank) + "----------")
print("SAVING COMPLETED IN:")
print(tsave)
#print("----------" + str(processorrank) + "----------")
