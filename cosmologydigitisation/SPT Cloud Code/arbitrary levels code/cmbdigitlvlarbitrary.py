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

# Read in data
filelocation = "../../for_lennart/plik_plus_r0p01_highell_lensedtotCls.dat"
file = open(filelocation)
l = []
dl = []
for row in file:
    rowvalues = [float(i) for i in row.split()]
    l.append(rowvalues[0])
    dl.append(rowvalues[1])

dl = asarray(dl)
l = asarray(l)


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

print(float(pixelsizearcmin)/raspeed)
fgsfes


# noisemap = 3.0 # muK arcmin
# noisepix = noisemap/float(pixelsizearcmin)
# fsky = (mapfieldsize*mapfieldsize)/(4.0*pi*60.0*60.0*(180.0/pi)**2.0)
# noisecl = 4.0*pi*fsky*noisepix*noisepix/float(mappixelnumber*mappixelnumber)
# eta = sqrt(float(compression)) * sqrt(float(observationlim)) * sqrt(float(noisecl) / (float(pixelnumber * pixelnumber)))
# eta0 = sqrt(float(observationlim)) * sqrt(float(noisecl) / (float(pixelnumber * pixelnumber)))


# Calculate cl
cl = (dl*2.0*pi)/(l*(l+1.0))
# Create 2dCl
cl2d = zeros((pixelnumber, pixelnumber))
pixelspacingrad = float(pixelsizearcmin)*arcmins2radians
kx = 2.0*pi*fftfreq(pixelnumber, pixelspacingrad)
kx2d = tile(kx, (pixelnumber, 1))
ky = 2.0*pi*fftfreq(pixelnumber, pixelspacingrad)
ky2d = transpose(tile(ky, (pixelnumber, 1)))
cl2d = zeros((pixelnumber, pixelnumber))
k2d = sqrt(kx2d*kx2d + ky2d*ky2d)
for y in range(pixelnumber):
    for x in range(pixelnumber):
        k = k2d[y, x]
        ind = argmin(abs(l-k))
        cl2d[y, x] = cl[ind]
#
# cl2d = zeros((pixelnumber, pixelnumber))
# fname = "../../cl2df" + str(fieldsizearcmins) + "r" + str(pixelsizearcmin) + ".txt"
# y = 0
# for row in open(fname):
#     rowvalues = row.split()
#     for x in range(pixelnumber):
#         cl2d[y, x] = float(rowvalues[x])
#     y += 1

#---

#df = 1.0/(float(fieldsizearcmins)*arcmins2radians)**2.0
cmbfft = normal(0, 1, (pixelnumber, pixelnumber)) + 1.0j * normal(0, 1, (pixelnumber, pixelnumber))
cmbfft = cmbfft * sqrt((df/2.0)*cl2d)
cmbmap = real(ifft2(cmbfft))[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
cmbmap = sqrt(2.0) * cmbmap * pixelnumber / (pixelsizearcmin * arcmins2radians)# * pixelnumber

#pixelrms = std(cmbmap)
#print(pixelrms)
#noisecl = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixelrms*pixelrms # muK^2
#print(noisecl)

#k, p, err, h = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin], cmbmap)
#p = p*pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians

imshow(cmbmap)
title(r"White noise $[\mu K]$", fontsize = 20)
colorbar()
show()

#plot(k, p, "b")
#plot(k, [noisecl for i in k], "r")
#title("Powerspectrum")
#xlabel(r"Multipole Moment $l$")
#ylabel(r"$C_l$")
#yscale("log")
#show()

#--- Define Noise

noisedet = 500.0 # muK sqrt(s)
noiseinduced = noisedet/(sqrt(1.0/float(readoutfreq))) # muK
pixelrms = noisedet/sqrt(float(pixelsizearcmin)/raspeed) # muK
noisemap = pixelrms*float(pixelsizearcmin) # muK arcmin
noisecl = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixelrms*pixelrms # muK^2
print(pixelrms)
print(noisecl)


#--- Add Noise




# #--- Read in CMP
# cmbmap = zeros((mappixelnumber, mappixelnumber))
# y = 0
# for row in open("cmbmaparblvltemplate.txt"):
#     rowvalues = row.split()
#     for x in range(mappixelnumber):
#         cmbmap[y, x] = float(rowvalues[x])
#     y += 1



#--- Define noise level

# noisedet = 500.0 # muK sqrt(s)
# noiseinduced = noisedet/(sqrt(1.0/float(readoutfreq))) # muK
# pixelrms = noisedet/sqrt(float(pixelsizearcmin)/raspeed) # muK
# noisemap = pixelrms*float(pixelsizearcmin) # muK arcmin
# #fsky = (mapfieldsize*mapfieldsize)/(4.0*pi*60.0*60.0*(180.0/pi)**2.0)
# #noisecl = 4.0*pi*fsky*pixelrms*pixelrms/float(mappixelnumber*mappixelnumber)
# noisecl = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixelrms*pixelrms # muK^2
#
#
# noisemap = normal(0, 2.0, (mappixelnumber, mappixelnumber))
#
#





# imshow(cmbmap)
# title("CMB MAP")
# colorbar()
# show()
#
# k, p, err, h = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin], cmbmap)
# dl = p*k*(k+1.0)/(pixelnumber*pixelnumber*pi)
# plot(k, dl)
# xlim((0, 2500))
# show()

#--- Recreate Observation
cmbnoisemap = zeros((nodecscans, norablocks))
#cmbnoisemap1bit = zeros((nodecscans, norablocks))
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
            #tod1bit = digitise1bit(tod, 1.0)
            #cmbnoisemap1bit[d, ri] += mean(tod1bit)


            #tod2bitopt = deepcopy(tod)

            #tod2bithm = digitise2bithalfmax(tod)
            #cmbnoisemap2bithm[d, ri] += mean(tod2bithm)

            #tod2bitopt = digitise2bitoptimal(tod2bitopt, eta)
            #cmbnoisemap2bitopt[d, ri] += mean(tod2bitopt)


        print("Noise, Digitisation and Compression: " + str(int(100*d/nodecscans)) + "%")

    print("--- COMPLETED OBSERVATION " + str(observationind) + " ---")

    observationind += 1

#--- SHOW MAPS

# subplot(1, 2, 1)
# imshow(cmbmap)
# title("CMB MAP")
# colorbar()
#
# subplot(1, 2, 2)
# imshow(cmbnoisemap1bit)
# title("CMB DIGIT MAP")
# colorbar()
# show()

#--- Normalise Maps

cmbnoisemap = cmbnoisemap * (1.0/float(observationlim))
#cmbnoisemap1bit = cmbnoisemap1bit * (1.0/float(observationlim))
#cmbnoisemap2bithm = cmbnoisemap2bithm * (1.0/float(observationlim))
#cmbnoisemap2bitopt = cmbnoisemap2bitopt * (1.0/float(observationlim))

#--- PS

k, p, err, h = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin], cmbnoisemap)
p = p * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians

plot(k, p, "b", label = "Reconstructed")
plot(l, dl*2.0*pi/(l*(l+1.0)), "r", label = "Input")
plot(k, [noisecl for o in k], "k", label = "Noise level")
xscale("log")
yscale("log")
title("CMB Powerspectrum", fontsize = 20)
xlabel(r"Multipole Moment $l$", fontsize = 20)
ylabel(r"$C_l [\mu K^2]$", fontsize = 20)
legend(loc = "lower left", fontsize = 15)
#xlim((0, 2500))
show()




# cmbmapadj = cmbmap*sqrt(sum(cmbnoisemap1bit**2.0)/sum(cmbmap**2.0))
# diffmap = cmbmapadj-cmbnoisemap1bit
#
# theo = linspace(-4, 4, 100)
# plot(cmbmapadj.flatten(), diffmap.flatten(), "kx")
# plot(theo, theo-(norm.cdf(theo)-norm.cdf(-theo)), "y")
# show()

#--- Save results


#savetxt("cmbmap" + str(processorrank) + ".txt", cmbmap)
#savetxt("cmbnoisemap" + str(processorrank+10) + ".txt", cmbnoisemap)
#savetxt("cmbnoisemap1bit" + str(processorrank) + ".txt", cmbnoisemap1bit)
#savetxt("cmbnoisemap2bithm" + str(processorrank+10) + ".txt", cmbnoisemap2bithm)
#savetxt("cmbnoisemap2bitopt" + str(processorrank+10) + ".txt", cmbnoisemap2bitopt)

#
# k0, p0, err0, h0 = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin], cmbmap)
k, p, err, h = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap)
# k1bit, p1bit, err1bit, h1bit = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap1bit)
# k2bithm, p2bithm, err2bithm, h2bithm = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap2bithm)
# k2bitopt, p2bitopt, err2bitopt, h2bitopt = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap2bitopt)

p = p/(pi * float(pixelnumber * pixelnumber))

p = p * k * (k + 1.0)

plot(linput, dlinput, "k")
plot(k, p, "b")
plot(k, [noisecl for i in k], "r")
#xscale("log")
yscale("log")
xlim((2500, max(k)))
show()




# # Normalise Powerspectrum
# mmi = 0
# truepow = 0
# obspow = 0
# obspow1bit = 0
# obspow2bithm = 0
# obspow2bitopt = 0
# while mmi < len(k0):
#     if k0[mmi] >= 350 and k0[mmi] <= 1000:
#         truepow += p0[mmi]**2.0
#         obspow += p[mmi]**2.0
#         obspow1bit += p1bit[mmi]**2.0
#         obspow2bithm += p2bithm[mmi]**2.0
#         obspow2bitopt += p2bitopt[mmi]**2.0
#     mmi += 1
#
# p = p * (truepow/obspow)**0.5
# p1bit = p1bit * (truepow/obspow1bit)**0.5
# p2bithm = p2bithm * (truepow/obspow2bithm)**0.5
# p2bitopt = p2bitopt * (truepow/obspow2bitopt)**0.5
#
# dltot = float(pixelnumber) * float(pixelnumber) * p0 * k0 * (k0 + 1.0) / pi
# dlobstot = float(pixelnumber) * float(pixelnumber) * p * k * (k + 1.0) / pi
# dlobstot1bit = float(pixelnumber) * float(pixelnumber) * p1bit * k1bit * (k1bit + 1.0) / pi
# dlobstot2bithm = float(pixelnumber) * float(pixelnumber) * p2bithm * k2bithm * (k2bithm + 1.0) / pi
# dlobstot2bitopt = float(pixelnumber) * float(pixelnumber) * p2bitopt * k2bitopt * (k2bitopt + 1.0) / pi
#
# cltot = (dltot*pi)/(k0*(k0+1.0))
# clobstot = (dlobstot*pi)/(k*(k+1.0))
# clobstot1bit = (dlobstot1bit*pi)/(k1bit*(k1bit+1.0))
# clobstot2bithm = (dlobstot2bithm*pi)/(k2bithm*(k2bithm+1.0))
# clobstot2bitopt = (dlobstot2bitopt*pi)/(k2bitopt*(k2bitopt+1.0))
#
# print("Saving results")
#
# savetxt("k0.txt", k0)
# savetxt("cltot.txt", cltot)
# savetxt("k.txt", k)
# savetxt("clobstot.txt", clobstot)
# savetxt("k1bit.txt", k1bit)
# savetxt("clobstot1bit.txt", clobstot1bit)
# savetxt("k2bithm.txt", k2bithm)
# savetxt("clobstot2bithm.txt", clobstot2bithm)
# savetxt("k2bitopt.txt", k2bitopt)
# savetxt("clobstot2bitopt.txt", clobstot2bitopt)
#
# print("Completed Calculation")