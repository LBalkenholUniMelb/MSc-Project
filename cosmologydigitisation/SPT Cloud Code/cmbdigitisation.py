#--- Necessary Imports
from numpy import zeros, shape, mean
from numpy.random import normal
from matplotlib.pyplot import subplot, imshow, title, colorbar, show

#--- Import and Define CMB
pixelnumber = 512
file = open("CMBMap.txt")
cmbmap = zeros((pixelnumber, pixelnumber))
progress = 0
rowindex = 0
for row in file:
    rowvalues = [float(i) for i in row.split()]
    cmbmap[rowindex] = rowvalues[:pixelnumber]
    progressnew = int(100 * rowindex / pixelnumber)
    if progress != progressnew:
        progress = progressnew
        print("Reading in CMB: " + str(int(100 * rowindex / pixelnumber)) + "%")
    rowindex += 1


var = 6*10**(-6) #K
avg = 0 #K

#--- Recreate Scan Strategy
declims = [0, 1024] #arcmins
ralims = [0, 1024] #arcmins
readoutfreq = 6 #Hz
raspeed = 0.1 #arcmin/s
nodecscans = 512
norablocks = 512
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
observationno = range(10)
observations = [zeros((nodecscans, norablocks)) for i in observationno]
cesscans = [zeros((nodecscans, radatapoints)) for i in observationno]

#--- Decompose into bolometer signals
noise = normal(avg, var, size = (len(observationno), nodecscans, norablocks*compression))
progress = 0
for d in range(nodecscans):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        # Create and noise
        for obs in observationno:
            tod = cmbmap[d, ri] + noise[obs][d][rstart:rstop]
            # digitise here
            #todpow = sum(asarray(tod)*asarray(tod))
            #digitise1bit(tod, 1)
            #toddigitpow = sum(asarray(tod)*asarray(tod))
            #tod = ((todpow/toddigitpow)**0.5)*tod
            cesscans[obs][d, rstart:rstop] = tod

    progressnew = int(100*d/nodecscans)
    if progress != progressnew:
        progress = progressnew
        print("Noise & Digitisation: " + str(int(100*d/nodecscans)) + "%")


#--- Recompress into map
progress = 0
for d in range(shape(cesscans[0])[0]):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        for obs in observationno:
            observations[obs][d, ri] = mean(cesscans[obs][d, rstart:rstop])

    progressnew = int(100 * d / nodecscans)
    if progress != progressnew:
        progress = progressnew
        print("Recompressing: " + str(int(100 * d / nodecscans)) + "%")

#--- Average Maps
cmbnoisemap = zeros((nodecscans, norablocks))
for obs in observations:
    cmbnoisemap = cmbnoisemap + obs
cmbnoisemap = cmbnoisemap * (1.0/len(observationno))

#--- Calculate further Statistics

#--- Plot/Output desired quantities
subplot(2, 2, 1)
imshow(cmbmap)
title("CMB")
colorbar()
subplot(2, 2, 2)
imshow(cmbnoisemap)
title("CMB + Noise")
colorbar()
subplot(2, 2, 3)
imshow(cmbmap-cmbnoisemap)
title("Difference")
colorbar()
show()