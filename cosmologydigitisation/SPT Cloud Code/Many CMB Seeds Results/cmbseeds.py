#--- Setup MPI
from mpi4py import MPI
comm = MPI.COMM_WORLD
processorrank = comm.Get_rank()

#--- Necessary Imports
from numpy import zeros, shape, mean, asarray, savetxt, arange, inf, sqrt, real
from numpy.fft import ifft2
from numpy.random import normal
from matplotlib.pyplot import subplot, imshow, title, colorbar, show, plot
from digitisationschemes import *


#--- Create CMB

# Define map parameters
fieldsizearcmins = 2048
pixelsizearcmin = 2
pixelnumber = 1024
df = 1.0

# Read in 2dCl
cl2d = zeros((pixelnumber, pixelnumber))
y = 0
for row in open("cl2d.txt"):
    rowvalues = row.split()
    for x in range(pixelnumber):
        cl2d[y, x] = float(rowvalues[x])
    y += 1


# Combine noise and cmb element-wise
factor = sqrt((df/2.0)*cl2d)

####################################################
####################################################
#------------ START LOOP FOR 10 CMB SEEDS HERE
####################################################
####################################################

for cmbseedindex in range(10):

    realpart = factor * normal(size = (pixelnumber, pixelnumber))
    imagpart = factor * normal(size = (pixelnumber, pixelnumber))
    cmbfreqspace = (realpart + 1j*imagpart)

    # Transform into map
    cmbmap = ifft2(cmbfreqspace)[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
    cmbmap = real(cmbmap)

    #--- Recreate Scan Strategy
    var = 6*10**(-6) #K
    avg = 0 #K
    declims = [0, 1024] #arcmins
    ralims = [0, 1024] #arcmins
    readoutfreq = 6.0 #Hz
    raspeed = 0.1 #arcmin/s
    nodecscans = 512
    norablocks = 512
    radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
    compression = int(radatapoints/norablocks)

    #--- Investigate multiple observation numbers
    allnobs = [0, 20]
    noise = normal(avg, var, size = (max(allnobs), nodecscans, norablocks*compression))

    allobservations1bit = [zeros((nodecscans, norablocks)) for i in range(max(allnobs))]
    allobservations2bithm = [zeros((nodecscans, norablocks)) for i in range(max(allnobs))]
    allobservations2bitopt = [zeros((nodecscans, norablocks)) for i in range(max(allnobs))]

    i = 0
    while i < len(allnobs)-1:
        formernobs = allnobs[i]
        nobs = allnobs[i+1]

        observations1bit = allobservations1bit[formernobs:nobs]
        observations2bithm = allobservations2bithm[formernobs:nobs]
        observations2bitopt = allobservations2bitopt[formernobs:nobs]

        observationno = range(nobs-formernobs)

        cesscans1bit = [zeros((nodecscans, radatapoints)) for j in observationno]
        cesscans2bithm = [zeros((nodecscans, radatapoints)) for j in observationno]
        cesscans2bitopt = [zeros((nodecscans, radatapoints)) for j in observationno]

        #--- Decompose into bolometer signals
        for d in range(nodecscans):
            for ri in range(norablocks):
                rstart = ri*compression
                rstop = rstart + compression
                # Create and add noise
                for obs in observationno:
                    tod1bit = cmbmap[d, ri] + noise[obs][d][rstart:rstop]
                    tod2bithm = cmbmap[d, ri] + noise[obs][d][rstart:rstop]
                    tod2bitopt = cmbmap[d, ri] + noise[obs][d][rstart:rstop]

                    todpow = sum(asarray(tod1bit)*asarray(tod1bit))

                    digitise1bit(tod1bit)
                    tod1bit = ((todpow/sum(tod1bit*tod1bit))** 0.5 ) * tod1bit
                    cesscans1bit[obs][d, rstart:rstop] = tod1bit

                    digitise2bithalfmax(tod2bithm)
                    tod2bithm = ((todpow / sum(tod2bithm * tod2bithm)) ** 0.5) * tod2bithm
                    cesscans2bithm[obs][d, rstart:rstop] = tod2bithm

                    digitise2bitoptimal(tod2bitopt, var)
                    tod2bitopt = ((todpow / sum(tod2bitopt * tod2bitopt)) ** 0.5) * tod2bitopt
                    cesscans2bitopt[obs][d, rstart:rstop] = tod2bitopt



        #--- Recompress into map
        for d in range(shape(cesscans1bit[0])[0]):
            for ri in range(norablocks):
                rstart = ri*compression
                rstop = rstart + compression
                for obs in observationno:

                    observations1bit[obs][d, ri] = mean(cesscans1bit[obs][d, rstart:rstop])
                    observations2bithm[obs][d, ri] = mean(cesscans2bithm[obs][d, rstart:rstop])
                    observations2bitopt[obs][d, ri] = mean(cesscans2bitopt[obs][d, rstart:rstop])

        #--- Insert recompressed map into observations list
        allobservations1bit[formernobs:nobs] = observations1bit
        allobservations2bithm[formernobs:nobs] = observations2bithm
        allobservations2bitopt[formernobs:nobs] = observations2bitopt

        i += 1


    #--- Average Maps and create desired output
    i = 0
    while i < len(allnobs)-1:
        formernobs = allnobs[i]
        nobs = allnobs[i+1]

        observations1bit = allobservations1bit[formernobs:nobs]
        observations2bithm = allobservations2bithm[formernobs:nobs]
        observations2bitopt = allobservations2bitopt[formernobs:nobs]

        cmbnoisemap1bit = zeros((nodecscans, norablocks))
        cmbnoisemap2bithm = zeros((nodecscans, norablocks))
        cmbnoisemap2bitopt = zeros((nodecscans, norablocks))

        f = 0
        while f < len(observations1bit):

            cmbnoisemap1bit = cmbnoisemap1bit + observations1bit[f]
            cmbnoisemap2bithm = cmbnoisemap2bithm + observations2bithm[f]
            cmbnoisemap2bitopt = cmbnoisemap2bitopt + observations2bitopt[f]

            f += 1

        cmbnoisemap1bit = cmbnoisemap1bit * 1.0/float(nobs)
        cmbnoisemap2bithm = cmbnoisemap2bithm * 1.0/float(nobs)
        cmbnoisemap2bitopt = cmbnoisemap2bitopt * 1.0/float(nobs)

        filename = "cmb" + "rank" + str(processorrank+12) + "seed" + str(cmbseedindex) + ".txt"
        savetxt(filename, cmbmap)
        filename = "digitisedcmb1bitnobs" + str(nobs) + "rank" + str(processorrank+12) + "seed" + str(cmbseedindex) + ".txt"
        savetxt(filename, cmbnoisemap1bit)
        filename = "digitisedcmb2bithmnobs" + str(nobs) + "rank" + str(processorrank+12) + "seed" + str(cmbseedindex) + ".txt"
        savetxt(filename, cmbnoisemap2bithm)
        filename = "digitisedcmb2bitnobsopt" + str(nobs) + "rank" + str(processorrank+12) + "seed" + str(cmbseedindex) + ".txt"
        savetxt(filename, cmbnoisemap2bitopt)

        i += 1

    ####################################################
    ####################################################
    #------------ FINISH LOOP FOR 10 CMB SEEDS HERE
    ####################################################
    ####################################################
