#--- Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from cosmdigitclasses import *
from numpy import *
from scipy.signal import convolve2d
from scipy.stats import *
from pickle import *


rc('text', usetex=True)
rc("xtick", labelsize = 15)
rc("ytick", labelsize = 15)

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
raspeed = 0.5 #arcmin/s
nodecscans = mappixelnumber
norablocks = mappixelnumber
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
observationlim = 10
compressions = [800000, 800000, 80000, 8000, 800, 80, 8]
raspeeds = [0.0005, 0.0005, 0.005, 0.05, 0.5, 5.0, 50.0]
obslims = [128, 10, 10, 10, 10, 10, 10]
mapparams = [mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin]
hpp = asarray(compressions)*asarray(obslims)
arcmins2radians = radians(1./60.)

#--- Define Noise

noisedet = 500.0 # muK sqrt(s)
noiseinduced = noisedet/(sqrt(1.0/float(readoutfreq))) # muK
pixelrms = noisedet/sqrt(float(pixelsizearcmin)/raspeed) # muK
pixelrmscombined = pixelrms/sqrt(float(observationlim))
noisemap = pixelrms*float(pixelsizearcmin) # muK arcmin
noisemapcombined = pixelrmscombined*float(pixelsizearcmin) # muK arcmin
noisecl = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixelrms*pixelrms # muK^2
noiseclcombined = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixelrmscombined*pixelrmscombined # muK^2

#--- READ IN MAPS

cmbmap = zeros((mappixelnumber, mappixelnumber))
y = 0
for row in open("cmbtemplatemapfixeddetnoise.txt"):
    rowvalues = row.split()
    for x in range(mappixelnumber):
        cmbmap[y, x] = float(rowvalues[x])
    y += 1

k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
p0 = p0 * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
dump(k0, open("k.p", "wb"), -1)
dump(p0, open("p0.p", "wb"), -1)

cmbmaps1bit = []
cmbmaps2bit = []
cmbmaps3bit = []
cmbnmaps1bit = []
cmbnmaps2bit = []
cmbnmaps3bit = []

i = 0
while i < len(hpp):

    if hpp[i] > 800000*100:
        cmbmaps1bit.append(load(open("results/cmbnoisemap1bitcombined.p", "rb")))
        cmbmaps2bit.append(load(open("results/cmbnoisemap2bitcombined.p", "rb")))
        cmbmaps3bit.append(load(open("results/cmbnoisemap3bitcombined.p", "rb")))
        cmbnmap = load(open("results/cmbnoisemapcombined.p", "rb"))
        cmbnmaps1bit.append(cmbnmap)
        cmbnmaps2bit.append(cmbnmap)
        cmbnmaps3bit.append(cmbnmap)

    else:

        cmbmap1bit = zeros((mappixelnumber, mappixelnumber))
        y = 0
        for row in open("../../SPT Cloud Code/fixed det noise results 1 bit/" + str(hpp[i]/10) + "hpp/cmbnoisemap1bitfixeddetnoise" + str(hpp[i]/10) + "hppcombined.txt", "rb"):
            rowvalues = row.split()
            for x in range(mappixelnumber):
                cmbmap1bit[y, x] = float(rowvalues[x])
            y += 1
        cmbmaps1bit.append(cmbmap1bit)

        cmbnmap1bit = zeros((mappixelnumber, mappixelnumber))
        y = 0
        for row in open("../../SPT Cloud Code/fixed det noise results 1 bit/" + str(hpp[i]/10) + "hpp/cmbnoisemapfixeddetnoise" + str(hpp[i]/10) + "hppcombined.txt", "rb"):
            rowvalues = row.split()
            for x in range(mappixelnumber):
                cmbnmap1bit[y, x] = float(rowvalues[x])
            y += 1
        cmbnmaps1bit.append(cmbnmap1bit)


        cmbmap2bit = zeros((mappixelnumber, mappixelnumber))
        y = 0
        for row in open("../../SPT Cloud Code/fixed det noise results 2 bit/" + str(hpp[i]/10) + "hpp/cmbnoisemap2bitoptfixeddetnoise" + str(hpp[i]/10) + "hppcombined.txt", "rb"):
            rowvalues = row.split()
            for x in range(mappixelnumber):
                cmbmap2bit[y, x] = float(rowvalues[x])
            y += 1
        cmbmaps2bit.append(cmbmap2bit)

        cmbnmap2bit = zeros((mappixelnumber, mappixelnumber))
        y = 0
        for row in open("../../SPT Cloud Code/fixed det noise results 2 bit/" + str(hpp[i]/10) + "hpp/cmbnoisemapfixeddetnoise" + str(hpp[i]/10) + "hppcombined.txt", "rb"):
            rowvalues = row.split()
            for x in range(mappixelnumber):
                cmbnmap2bit[y, x] = float(rowvalues[x])
            y += 1
        cmbnmaps2bit.append(cmbnmap2bit)


        cmbmaps3bit.append(load(open("../../SPT Cloud Code/fixed det noise results 3 bit/" + str(hpp[i]/10) + "hpp/cmbnoisemap3bitfixeddetnoise" + str(hpp[i]/10) + "hppcombined.p", "rb")))

        cmbnmaps3bit.append(load(open("../../SPT Cloud Code/fixed det noise results 3 bit/" + str(hpp[i]/10) + "hpp/cmbnoisemap" + str(hpp[i]/10) + "hppcombined.p", "rb")))


    i += 1

print("#--------------------#")
print("COMPLETED READING IN DATA")
print("#--------------------#")


#--- CALCULATE AND NORMALISE POWERSPECTRA


k = 0
ps1bit = []
psn1bit = []
ps2bit = []
psn2bit = []
ps3bit = []
psn3bit = []

###
k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
p0 = p0 * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
dump(k, open("k.p", "wb"), -1)
dump(p0, open("p0.p", "wb"), -1)

i = 0
while i < len(hpp):
    k1bit, p1bit, err1bit, h1bit = cosm.sriniPowerSpectrum(mapparams, cmbmaps1bit[i])
    k2bit, p2bit, err2bit, h2bit = cosm.sriniPowerSpectrum(mapparams, cmbmaps2bit[i])
    k3bit, p3bit, err3bit, h3bit = cosm.sriniPowerSpectrum(mapparams, cmbmaps3bit[i])
    kn1bit, pn1bit, errn1bit, hn1bit = cosm.sriniPowerSpectrum(mapparams, cmbnmaps1bit[i])
    kn2bit, pn2bit, errn2bit, hn2bit = cosm.sriniPowerSpectrum(mapparams, cmbnmaps2bit[i])
    kn3bit, pn3bit, errn3bit, hn3bit = cosm.sriniPowerSpectrum(mapparams, cmbnmaps3bit[i])
    p1bit = p1bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    p2bit = p2bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    p3bit = p3bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    pn1bit = pn1bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    pn2bit = pn2bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    pn3bit = pn3bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians

    ps1bit.append(p1bit)
    psn1bit.append(pn1bit)
    ps2bit.append(p2bit)
    psn2bit.append(pn2bit)
    ps3bit.append(p3bit)
    psn3bit.append(pn3bit)

    if i == 0:
        k = k1bit
    i += 1

print("#--------------------#")
print("COMPLETED CALCULATING POWERSPECTRA")
print("#--------------------#")


normwindows = [[350, 1000], [350, 1000], [350, 1000], [350, 1000], [350, 800], [50, 300], [15, 100]]

i = 0
while i < len(hpp):

    p1bit = ps1bit[i]
    pn1bit = psn1bit[i]
    p2bit = ps2bit[i]
    pn2bit = psn2bit[i]
    p3bit = ps3bit[i]
    pn3bit = psn3bit[i]

    mmi = 0
    obspow3 = 0
    obspow2 = 0
    obspow1 = 0
    obspow3bit = 0
    obspow2bit = 0
    obspow1bit = 0

    while mmi < len(k):
        if k[mmi] >= normwindows[i][0] and k[mmi] <= normwindows[i][1]:
            obspow3 += pn3bit[mmi] ** 2.0
            obspow2 += pn2bit[mmi] ** 2.0
            obspow1 += pn1bit[mmi] ** 2.0
            obspow3bit += p3bit[mmi] ** 2.0
            obspow2bit += p2bit[mmi] ** 2.0
            obspow1bit += p1bit[mmi] ** 2.0
        mmi += 1

    p3bit = p3bit * (obspow3 / obspow3bit) ** 0.5
    p2bit = p2bit * (obspow2 / obspow2bit) ** 0.5
    p1bit = p1bit * (obspow1 / obspow1bit) ** 0.5

    ps1bit[i] = p1bit
    ps2bit[i] = p2bit
    ps3bit[i] = p3bit

    i += 1

print("#--------------------#")
print("COMPLETED NORMALISING POWERSPECTRA")
print("#--------------------#")


dump(k, open("k.p", "wb"), -1)
dump(p0, open("p0.p", "wb"), -1)

i = 0
while i < len(hpp):

    p1bit = ps1bit[i]
    pn1bit = psn1bit[i]
    p2bit = ps2bit[i]
    pn2bit = psn2bit[i]
    p3bit = ps3bit[i]
    pn3bit = psn3bit[i]

    dump(p1bit, open("p1bit" + str(int(hpp[i])) + "hpp.p", "wb"), -1)
    dump(pn1bit, open("pn1bit" + str(int(hpp[i])) + "hpp.p", "wb"), -1)
    dump(p2bit, open("p2bit" + str(int(hpp[i])) + "hpp.p", "wb"), -1)
    dump(pn2bit, open("pn2bit" + str(int(hpp[i])) + "hpp.p", "wb"), -1)
    dump(p3bit, open("p3bit" + str(int(hpp[i])) + "hpp.p", "wb"), -1)
    dump(pn3bit, open("pn3bit" + str(int(hpp[i])) + "hpp.p", "wb"), -1)

    i += 1