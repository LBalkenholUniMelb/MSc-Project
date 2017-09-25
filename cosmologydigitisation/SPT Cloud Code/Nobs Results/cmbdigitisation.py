#--- Necessary Imports
from numpy import zeros, shape, mean, asarray, savetxt, arange, inf
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
allnobs = [0, 20]#[0, 1] + list(arange(10, 210, 10))
noise = normal(avg, var, size = (max(allnobs), nodecscans, norablocks*compression))

allobservations1bit = [zeros((nodecscans, norablocks)) for i in range(max(allnobs))]
allobservations2bithm = [zeros((nodecscans, norablocks)) for i in range(max(allnobs))]
allobservations2biteqn = [zeros((nodecscans, norablocks)) for i in range(max(allnobs))]
allobservations2bitopt = [zeros((nodecscans, norablocks)) for i in range(max(allnobs))]

i = 0
while i < len(allnobs)-1:
    formernobs = allnobs[i]
    nobs = allnobs[i+1]

    observations1bit = allobservations1bit[formernobs:nobs]
    observations2bithm = allobservations2bithm[formernobs:nobs]
    observations2biteqn = allobservations2biteqn[formernobs:nobs]
    observations2bitopt = allobservations2bitopt[formernobs:nobs]

    observationno = range(nobs-formernobs)

    cesscans1bit = [zeros((nodecscans, radatapoints)) for j in observationno]
    cesscans2bithm = [zeros((nodecscans, radatapoints)) for j in observationno]
    cesscans2biteqn = [zeros((nodecscans, radatapoints)) for j in observationno]
    cesscans2bitopt = [zeros((nodecscans, radatapoints)) for j in observationno]

    #--- Decompose into bolometer signals
    progress = 0
    for d in range(nodecscans):
        for ri in range(norablocks):
            rstart = ri*compression
            rstop = rstart + compression
            # Create and add noise
            for obs in observationno:
                tod1bit = cmbmap[d, ri] + noise[obs][d][rstart:rstop]
                tod2bithm = cmbmap[d, ri] + noise[obs][d][rstart:rstop]
                tod2biteqn = cmbmap[d, ri] + noise[obs][d][rstart:rstop]
                tod2bitopt = cmbmap[d, ri] + noise[obs][d][rstart:rstop]

                todpow = sum(asarray(tod1bit)*asarray(tod1bit))

                digitise1bit(tod1bit)
                tod1bit = ( (todpow/sum(tod1bit*tod1bit))**0.5 ) * tod1bit
                cesscans1bit[obs][d, rstart:rstop] = tod1bit

                digitise2bithalfmax(tod2bithm)
                tod2bithm = ((todpow / sum(tod2bithm * tod2bithm)) ** 0.5) * tod2bithm
                cesscans2bithm[obs][d, rstart:rstop] = tod2bithm

                digitise2bitequalnumbers(tod2biteqn)
                tod2biteqn = ((todpow / sum(tod2biteqn * tod2biteqn)) ** 0.5) * tod2biteqn
                cesscans2biteqn[obs][d, rstart:rstop] = tod2biteqn

                digitise2bitoptimal(tod2bitopt)
                tod2bitopt = ((todpow / sum(tod2bitopt * tod2bitopt)) ** 0.5) * tod2bitopt
                cesscans2bitopt[obs][d, rstart:rstop] = tod2bitopt

        progressnew = 10*int(10*d/nodecscans)
        if progress != progressnew:
            progress = progressnew
            print("Noise & Digitisation: " + str(10*int(10*d/nodecscans)) + "%")


    #--- Recompress into map
    progress = 0
    for d in range(shape(cesscans1bit[0])[0]):
        for ri in range(norablocks):
            rstart = ri*compression
            rstop = rstart + compression
            for obs in observationno:

                observations1bit[obs][d, ri] = mean(cesscans1bit[obs][d, rstart:rstop])
                observations2bithm[obs][d, ri] = mean(cesscans2bithm[obs][d, rstart:rstop])
                observations2biteqn[obs][d, ri] = mean(cesscans2biteqn[obs][d, rstart:rstop])
                observations2bitopt[obs][d, ri] = mean(cesscans2bitopt[obs][d, rstart:rstop])

        progressnew = 10*int(10 * d / nodecscans)
        if progress != progressnew:
            progress = progressnew
            print("Recompressing: " + str(10*int(10 * d / nodecscans)) + "%")

    #--- Insert recompressed map into observations list
    allobservations1bit[formernobs:nobs] = observations1bit
    allobservations2bithm[formernobs:nobs] = observations2bithm
    allobservations2biteqn[formernobs:nobs] = observations2biteqn
    allobservations2bitopt[formernobs:nobs] = observations2bitopt

    print("------------------------------")
    print("Completed: " + str(nobs) + " Nobs")
    print("------------------------------")
    i += 1


#--- Average Maps and create desired output
i = 0
while i < len(allnobs)-1:
    formernobs = allnobs[i]
    nobs = allnobs[i+1]

    observations1bit = allobservations1bit[formernobs:nobs]
    observations2bithm = allobservations2bithm[formernobs:nobs]
    observations2biteqn = allobservations2biteqn[formernobs:nobs]
    observations2bitopt = allobservations2bitopt[formernobs:nobs]

    cmbnoisemap1bit = zeros((nodecscans, norablocks))
    cmbnoisemap2bithm = zeros((nodecscans, norablocks))
    cmbnoisemap2biteqn = zeros((nodecscans, norablocks))
    cmbnoisemap2bitopt = zeros((nodecscans, norablocks))

    f = 0
    while f < len(observations1bit):
        cmbnoisemap1bit = cmbnoisemap1bit + observations1bit[f]
        cmbnoisemap2bithm = cmbnoisemap2bithm + observations2bithm[f]
        cmbnoisemap2biteqn = cmbnoisemap2biteqn + observations2biteqn[f]
        cmbnoisemap2bitopt = cmbnoisemap2bitopt + observations2bitopt[f]
        f += 1

    cmbnoisemap1bit = cmbnoisemap1bit * float(1/nobs)
    cmbnoisemap2bithm = cmbnoisemap2bithm * float(1/nobs)
    cmbnoisemap2biteqn = cmbnoisemap2biteqn * float(1/nobs)
    cmbnoisemap2bitopt = cmbnoisemap2bitopt * float(1/nobs)

    filename = "digitisedcmb1bitnobs" + str(nobs) + ".txt"
    savetxt(filename, cmbnoisemap1bit)
    filename = "digitisedcmb2bithmnobs" + str(nobs) + ".txt"
    savetxt(filename, cmbnoisemap2bithm)
    filename = "digitisedcmb2biteqnnobs" + str(nobs) + ".txt"
    savetxt(filename, cmbnoisemap2biteqn)
    filename = "digitisedcmb2bitnobsopt" + str(nobs) + ".txt"
    savetxt(filename, cmbnoisemap2bitopt)

    i += 1