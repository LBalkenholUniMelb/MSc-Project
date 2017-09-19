#--- Necessary Imports
from numpy import zeros, shape, mean, asarray, savetxt, arange
from numpy.random import normal
from matplotlib.pyplot import subplot, imshow, title, colorbar, show
from digitisationschemes import digitise1bit

#--- Import and Define CMB
pixelnumber = 512
file = open("CMBMap.txt")
cmbmap = zeros((pixelnumber, pixelnumber))
rowindex = 0
for row in file:
    rowvalues = [float(i) for i in row.split()]
    cmbmap[rowindex] = rowvalues[:pixelnumber]
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

#--- Investigate multiple observation numbers
allnobs = [0, 1] + list(arange(10, 210, 10))
noise = normal(avg, var, size = (max(allnobs), nodecscans, norablocks*compression))
allobservations = [zeros((nodecscans, norablocks)) for i in range(max(allnobs))]
i = 0
while i < len(allnobs)-1:
    formernobs = allnobs[i]
    nobs = allnobs[i+1]
    observations = allobservations[formernobs:nobs]
    observationno = range(nobs-formernobs)
    cesscans = [zeros((nodecscans, radatapoints)) for j in observationno]

    #--- Decompose into bolometer signals
    progress = 0
    for d in range(nodecscans):
        for ri in range(norablocks):
            rstart = ri*compression
            rstop = rstart + compression
            # Create and add noise
            for obs in observationno:
                tod = cmbmap[d, ri] + noise[obs][d][rstart:rstop]
                # digitise here
                todpow = sum(asarray(tod)*asarray(tod))
                digitise1bit(tod)
                toddigitpow = sum(asarray(tod)*asarray(tod))
                tod = ((todpow/toddigitpow)**0.5)*tod
                cesscans[obs][d, rstart:rstop] = tod

        progressnew = 10*int(10*d/nodecscans)
        if progress != progressnew:
            progress = progressnew
            print("Noise & Digitisation: " + str(10*int(10*d/nodecscans)) + "%")


    #--- Recompress into map
    progress = 0
    for d in range(shape(cesscans[0])[0]):
        for ri in range(norablocks):
            rstart = ri*compression
            rstop = rstart + compression
            for obs in observationno:
                observations[obs][d, ri] = mean(cesscans[obs][d, rstart:rstop])

        progressnew = 10*int(10 * d / nodecscans)
        if progress != progressnew:
            progress = progressnew
            print("Recompressing: " + str(10*int(10 * d / nodecscans)) + "%")

    #--- Insert recompressed map into observations list
    allobservations[formernobs:nobs] = observations
    print("------------------------------")
    print("Completed: " + str(nobs) + " Nobs")
    print("------------------------------")
    i += 1

#--- Average Maps and create desired output
i = 0
while i < len(allnobs)-1:
    formernobs = allnobs[i]
    nobs = allnobs[i+1]
    observations = allobservations[formernobs:nobs]
    cmbnoisemap = zeros((nodecscans, norablocks))
    for obs in observations:
        cmbnoisemap = cmbnoisemap + obs
    cmbnoisemap = cmbnoisemap * (1.0/float(nobs))
    filename = "digitisedcmbnobs" + str(nobs) + ".txt"
    savetxt(filename, cmbnoisemap)
    i += 1