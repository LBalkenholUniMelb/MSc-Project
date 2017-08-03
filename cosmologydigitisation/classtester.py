from cosmdigitclasses import *
from numpy import corrcoef, average, savetxt
#from pow_spec import *

# Observation Mock Data
observationstarttime = 0
observationendtime = 10
readoutfrequency = 200
numberofreadouts = readoutfrequency*(observationendtime-observationstarttime)
numberofobservations = 200

# Signal Mock Data
var = 0.0006
meanvalue = 2.7255

# Map Mock Data
raresolution = 100
decresolution = 100
philim = [0.0, 2.0]
thetalim = [0.0, 2.0]

# Create Mock Data

#signaldigit = Digitiser()
mockobservationspure = []
#mockobservations1bit = []
#mockobservations2bitopt = []
#mockobservations2bithalf = []
#mockobservations2biteq = []

# Compare noise signals to digitised signals
#gwn = GaussianWhiteNoise(observationendtime - observationstarttime, readoutfrequency, meanvalue, var)
cosm = Cosmologist()
# cosm.compareGwn1Bit(gwn, colours=[["k"], ["r"]])
# cosm.compareGwn2BitDigit(gwn, scheme="halfmax", colours=[["k"], ["r"]])
# cosm.compareGwn2BitDigit(gwn, scheme="equalnumbers", colours=[["k"], ["r"]])
# cosm.compareGwn2BitDigit(gwn, scheme="optimal", colours=[["k"], ["r"]])

for i in range(numberofobservations):
    # Create Noise Signal and Digitise According to Different Schemes
    gwn = GaussianWhiteNoise(observationendtime - observationstarttime, readoutfrequency, meanvalue, var)
    #gwn1bit = signaldigit.digitise1bit(gwn)

#     gwn1bit = signaldigit.digitise1bit(gwn, meanvalue)
#     gwn2bitopt = signaldigit.digitise2bit(gwn, "optimal", meanvalue, var)
#     gwn2bithalf = signaldigit.digitise2bit(gwn, "halfmax", meanvalue, var)
#     gwn2biteq = signaldigit.digitise2bit(gwn, "equalnumbers", meanvalue, var)
#
    # Create Observations
    phiscan = linspace(philim[0], philim[1], numberofreadouts)
    thetastep = (thetalim[1]-thetalim[0])/numberofobservations
    theta = [(thetalim[0] + (thetastep/2) + (thetastep*i)) for k in range(numberofreadouts)]
    mockobservationspure.append(Observation(observationstarttime, observationendtime, readoutfrequency, phiscan, theta, gwn))
    #mockobservations1bit.append(Observation(observationstarttime, observationendtime, readoutfrequency, phiscan, theta, gwn1bit))

#     mockobservations1bit.append(Observation(observationstarttime, observationendtime, readoutfrequency, phiscan, theta, gwn1bit))
#     mockobservations2bitopt.append(Observation(observationstarttime, observationendtime, readoutfrequency, phiscan, theta, gwn2bitopt))
#     mockobservations2bithalf.append(Observation(observationstarttime, observationendtime, readoutfrequency, phiscan, theta, gwn2bithalf))
#     mockobservations2biteq.append(Observation(observationstarttime, observationendtime, readoutfrequency, phiscan, theta, gwn2bithalf))

# Create Maps from Mock Data
mockmappure = Map(mockobservationspure, raresolution, decresolution)
mockmappure.plotMap(noise = True)
#mockmap1bit = Map(mockobservations1bit, raresolution, decresolution)
#mockmap1bit.plotMap(noise = True)

cosm.sriniPowerSpectrum([raresolution, decresolution, int(philim[1]), int(thetalim[1])], mockmappure)
#cosm.sriniPowerSpectrum([raresolution, decresolution, int(philim[1]), int(thetalim[1])], mockmap1bit)


#mockmap1bit = Map(mockobservations1bit, raresolution, decresolution)
#mockmap2bitopt = Map(mockobservations2bitopt, raresolution, decresolution)
#mockmap2bithalf = Map(mockobservations2bithalf, raresolution, decresolution)
#mockmap2biteq = Map(mockobservations2biteq, raresolution, decresolution)

# Calculate Powerspectra
# calculate map parameters
# nx, ny = [int(philim[1]-philim[0]), int(thetalim[1]-thetalim[0])]
# dx = (nx/raresolution)*60
# dy = (ny/decresolution)*60
# powspecpure = fn_plot_pow_spec([nx, ny, dx, dy], mockmappure.map)
# print(powspecpure)
# #powspec1bit = fn_plot_pow_spec([nx, ny, dx, dy], mockmap1bit.map)
# #powspec2bitopt = fn_plot_pow_spec([nx, ny, dx, dy], mockmap2bitopt.map)
# #powspec2bithalf = fn_plot_pow_spec([nx, ny, dx, dy], mockmap2bithalf.map)
# #powspec2biteq = fn_plot_pow_spec([nx, ny, dx, dy], mockmap2biteq.map)
#
# errorbar(powspecpure[0][1:], powspecpure[1][1:], yerr = powspecpure[2][1:], label = "Pure")
# #errorbar(powspec1bit[0][1:], powspec1bit[1][1:], yerr = powspec1bit[2][1:], label = "1 Bit")
# #errorbar(powspec2bitopt[0][1:], powspec2bitopt[1][1:], yerr = powspec2bitopt[2][1:], label = "2 Bit Opt")
# #errorbar(powspec2bithalf[0][1:], powspec2bithalf[1][1:], yerr = powspec2bithalf[2][1:], label = "2 Bit Half")
# #errorbar(powspec2biteq[0][1:], powspec2biteq[1][1:], yerr = powspec2biteq[2][1:], label = "2 Bit Eq")
#
# # plot cross power
# title("Power Spectra, GWN, digitised")
# xlabel("l")
# ylabel("pl")
# legend()
# show()