# Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from scipy.stats import binned_statistic
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from cosmdigitclasses import *
from numpy import *
from scipy.signal import convolve2d
rc('text', usetex=True)
rc("xtick", labelsize = 15)
rc("ytick", labelsize = 15)

cosm = Cosmologist()

# Define map parameters
fieldsizearcmins = 2048
pixelsizearcmin = 2
pixelnumber = int(fieldsizearcmins/pixelsizearcmin)
df = 1.0
mapfieldsize = int(fieldsizearcmins/2.0)
mappixelnumber = int(pixelnumber/2.0)
declims = [0, mapfieldsize] #arcmins
ralims = [0, mapfieldsize] #arcmins
readoutfreq = 40 #Hz
raspeed = 0.01 #arcmin/s
nodecscans = mappixelnumber
norablocks = mappixelnumber
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
observationlim = 30

print(compression)
fmseljfe

# noise level
noisemap = 3.0 # muK arcmin
noisepix = noisemap/float(pixelsizearcmin)
fsky = (mapfieldsize*mapfieldsize)/(4.0*pi*60.0*60.0*(180.0/pi)**2.0)
noisecl = 4.0*pi*fsky*noisepix*noisepix/float(mappixelnumber*mappixelnumber)
eta = sqrt(float(compression)) * sqrt(float(observationlim)) * sqrt(float(noisecl) / (float(pixelnumber * pixelnumber)))


# Read in 2dCl
cl2d = zeros((pixelnumber, pixelnumber))
fname = "cl2df" + str(fieldsizearcmins) + "r" + str(pixelsizearcmin) + ".txt"
y = 0
for row in open(fname):
    rowvalues = row.split()
    for x in range(pixelnumber):
        cl2d[y, x] = float(rowvalues[x])
    y += 1

# Combine noise and cmb element-wise

factor = sqrt((df/2.0)*cl2d)
facadj = 2.0

realno = 1
realind = 0
k0 = 0
k = 0
dltot = 0
dlobstot = 0
dlobstot1bit = 0
dlobstot2bithm = 0
dlobstot2bitopt = 0

while realind < realno:

    re = normal(0, 1, (pixelnumber, pixelnumber))
    im = normal(0, 1, (pixelnumber, pixelnumber))

    realpart = factor * re
    imagpart = factor * im
    cmbfreqspace = (realpart + 1.0j*imagpart)

    # Transform into map
    cmbmap = fft.ifft2(cmbfreqspace)[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
    cmbmap = real(cmbmap)

    # Recreate Scan Strategy


    observationno = range(observationlim)
    observations = [zeros((nodecscans, norablocks)) for i in observationno]
    observations1bit = [zeros((nodecscans, norablocks)) for i in observationno]
    observations2bithm = [zeros((nodecscans, norablocks)) for i in observationno]
    observations2bitopt = [zeros((nodecscans, norablocks)) for i in observationno]

    cesscans = [zeros((nodecscans, radatapoints)) for i in observationno]
    cesscans1bit = [zeros((nodecscans, radatapoints)) for i in observationno]
    cesscans2bithm = [zeros((nodecscans, radatapoints)) for i in observationno]
    cesscans2bitopt = [zeros((nodecscans, radatapoints)) for i in observationno]

    # Decompose into bolometer signals
    for d in range(nodecscans):

        # create noise for this scan
        noiserall = [0 for pp in range(observationlim)]
        noisevarall = [0 for uu in range(observationlim)]
        noiseind = 0
        while noiseind < observationlim:
            noisef = normal(0, eta ** 0.5, norablocks * compression)
            noisef = noisef
            noiser = ifft(noisef)
            noiser = real(noiser)
            noiserall[noiseind] = noiser
            noisevarall[noiseind] = std(noiser)
            noiseind += 1


        for ri in range(norablocks):
            rstart = ri*compression
            rstop = rstart + compression

            for obs in observationno:

                # add noise
                tod = cmbmap[d, ri] + noiserall[obs][rstart:rstop]
                tod = asarray([cmbmap[d, ri] for x in range(compression)])
                cesscans[obs][d, rstart:rstop] = tod

                # digitise here
                tod1bit = digitise1bit(tod, 1.0)
                cesscans1bit[obs][d, rstart:rstop] = tod1bit

                tod2bithm = digitise2bithalfmax(tod)
                cesscans2bithm[obs][d, rstart:rstop] = tod2bithm

                tod2bitopt = digitise2bitoptimal(tod, noisevarall[obs])
                cesscans2bitopt[obs][d, rstart:rstop] = tod2bitopt


        print("Noise & Digitisation: " + str(int(100*d/nodecscans)) + "%")

    # Recompress into map

    print("Decomposition completed")


    for d in range(shape(cesscans[0])[0]):
        for ri in range(norablocks):
            rstart = ri*compression
            rstop = rstart + compression
            for obs in observationno:
                observations[obs][d, ri] = mean(cesscans[obs][d, rstart:rstop])
                observations1bit[obs][d, ri] = mean(cesscans1bit[obs][d, rstart:rstop])
                observations2bithm[obs][d, ri] = mean(cesscans2bithm[obs][d, rstart:rstop])
                observations2bitopt[obs][d, ri] = mean(cesscans2bitopt[obs][d, rstart:rstop])


        print("Compression: " + str(int(100*d/nodecscans)) + "%")


    print("Beginning PS extraction")

    cmbnoisemap = zeros((nodecscans, norablocks))
    cmbnoisemap1bit = zeros((nodecscans, norablocks))
    cmbnoisemap2bithm = zeros((nodecscans, norablocks))
    cmbnoisemap2bitopt = zeros((nodecscans, norablocks))

    for obs in range(observationlim):
        cmbnoisemap = cmbnoisemap + observations[obs]
        cmbnoisemap1bit = cmbnoisemap1bit + observations1bit[obs]
        cmbnoisemap2bithm = cmbnoisemap2bithm + observations2bithm[obs]
        cmbnoisemap2bitopt = cmbnoisemap2bitopt + observations2bitopt[obs]


    cmbnoisemap = cmbnoisemap * (1.0/len(observationno))
    cmbnoisemap1bit = cmbnoisemap1bit * (1.0/len(observationno))
    cmbnoisemap2bithm = cmbnoisemap2bithm * (1.0/len(observationno))
    cmbnoisemap2bitopt = cmbnoisemap2bitopt * (1.0/len(observationno))


    k0, p0, err0, h0 = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin], cmbmap)
    k, p, err, h = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap)
    k1bit, p1bit, err1bit, h1bit = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap1bit)
    k2bithm, p2bithm, err2bithm, h2bithm = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap2bithm)
    k2bitopt, p2bitopt, err2bitopt, h2bitopt = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap2bitopt)


    # Normalise Powerspectrum
    mmi = 0
    truepow = 0
    obspow = 0
    obspow1bit = 0
    obspow2bithm = 0
    obspow2bitopt = 0
    while mmi < len(k0):
        if k0[mmi] >= 350 and k0[mmi] <= 1000:
            truepow += p0[mmi]**2.0
            obspow += p[mmi]**2.0
            obspow1bit += p1bit[mmi]**2.0
            obspow2bithm += p2bithm[mmi]**2.0
            obspow2bitopt += p2bitopt[mmi]**2.0
        mmi += 1

    p = p * (truepow/obspow)**0.5
    p1bit = p1bit * (truepow/obspow1bit)**0.5
    p2bithm = p2bithm * (truepow/obspow2bithm)**0.5
    p2bitopt = p2bitopt * (truepow/obspow2bitopt)**0.5

    dltot = float(pixelnumber) * float(pixelnumber) * p0 * k0 * (k0 + 1.0) / pi
    dlobstot = float(pixelnumber) * float(pixelnumber) * p * k * (k + 1.0) / pi
    dlobstot1bit = float(pixelnumber) * float(pixelnumber) * p1bit * k1bit * (k1bit + 1.0) / pi
    dlobstot2bithm = float(pixelnumber) * float(pixelnumber) * p2bithm * k2bithm * (k2bithm + 1.0) / pi
    dlobstot2bitopt = float(pixelnumber) * float(pixelnumber) * p2bitopt * k2bitopt * (k2bitopt + 1.0) / pi


    realind += 1



cltot = (dltot*pi)/(k0*(k0+1.0))
clobstot = (dlobstot*pi)/(k*(k+1.0))
clobstot1bit = (dlobstot1bit*pi)/(k1bit*(k1bit+1.0))
clobstot2bithm = (dlobstot2bithm*pi)/(k2bithm*(k2bithm+1.0))
clobstot2bitopt = (dlobstot2bitopt*pi)/(k2bitopt*(k2bitopt+1.0))

print("Saving results")

savetxt("k0.txt", k0)
savetxt("cltot.txt", cltot)
savetxt("k.txt", k)
savetxt("clobstot.txt", clobstot)
savetxt("k1bit.txt", k1bit)
savetxt("clobstot1bit.txt", clobstot1bit)
savetxt("k2bithm.txt", k2bithm)
savetxt("clobstot2bithm.txt", clobstot2bithm)
savetxt("k2bitopt.txt", k2bitopt)
savetxt("clobstot2bitopt.txt", clobstot2bitopt)

print("Completed Calculation")