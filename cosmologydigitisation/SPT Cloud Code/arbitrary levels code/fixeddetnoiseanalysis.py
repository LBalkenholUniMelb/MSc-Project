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
compressions = [800000, 80000, 8000, 800, 80, 8]
raspeeds = [0.0005, 0.005, 0.05, 0.5, 5.0, 50.0]


#--- Define Noise

noisedet = 500.0 # muK sqrt(s)
noiseinduced = noisedet/(sqrt(1.0/float(readoutfreq))) # muK
pixelrms = noisedet/sqrt(float(pixelsizearcmin)/raspeed) # muK
pixelrmscombined = pixelrms/sqrt(float(observationlim))
noisemap = pixelrms*float(pixelsizearcmin) # muK arcmin
noisemapcombined = pixelrmscombined*float(pixelsizearcmin) # muK arcmin
noisecl = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixelrms*pixelrms # muK^2
noiseclcombined = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixelrmscombined*pixelrmscombined # muK^2



# #-----------------------#
# #--- 3 Bit interlude ---#
# #-----------------------#
#
#
# mapparams = [mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin]
# cmbmap = zeros((mappixelnumber, mappixelnumber))
# y = 0
# for row in open("../../NERSC Code/fixeddetnoise/cmbtemplatemapfixeddetnoise.txt"):
#     rowvalues = row.split()
#     for x in range(mappixelnumber):
#         cmbmap[y, x] = float(rowvalues[x])
#     y += 1
# cmbnmap = load(open("../../NERSC Code/fixeddetnoise/results/cmbnoisemapcombined.p", "rb"))
# cmbnmap1bit = load(open("../../NERSC Code/fixeddetnoise/results/cmbnoisemap1bitcombined.p", "rb"))
# cmbnmap2bit = load(open("../../NERSC Code/fixeddetnoise/results/cmbnoisemap2bitcombined.p", "rb"))
# cmbnmap3bit = load(open("../../NERSC Code/fixeddetnoise/results/cmbnoisemap3bitcombined.p", "rb"))
#
# k, p, err, h = cosm.sriniPowerSpectrum(mapparams, cmbmap)
# kn, pn, errn, hn = cosm.sriniPowerSpectrum(mapparams, cmbnmap)
# kn1bit, pn1bit, errn1bit, hn1bit = cosm.sriniPowerSpectrum(mapparams, cmbnmap1bit)
# kn2bit, pn2bit, errn2bit, hn2bit = cosm.sriniPowerSpectrum(mapparams, cmbnmap2bit)
# kn3bit, pn3bit, errn3bit, hn3bit = cosm.sriniPowerSpectrum(mapparams, cmbnmap3bit)
#
# p = p * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
# pn = pn * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
# pn1bit = pn1bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
# pn2bit = pn2bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
# pn3bit = pn3bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
#
# mmi = 0
# obspow = 0
# obspow1bit = 0
# obspow2bit = 0
# obspow3bit = 0
# while mmi < len(kn):
#     if kn[mmi] >= 350 and kn[mmi] <= 1000:
#         obspow += pn[mmi] ** 2.0
#         obspow1bit += pn1bit[mmi] ** 2.0
#         obspow2bit += pn2bit[mmi] ** 2.0
#         obspow3bit += pn3bit[mmi] ** 2.0
#     mmi += 1
# pn1bit = pn1bit * (obspow / obspow1bit) ** 0.5
# pn2bit = pn2bit * (obspow / obspow2bit) ** 0.5
# pn3bit = pn3bit * (obspow / obspow3bit) ** 0.5
#
# indstart = argmin(abs(kn1bit - 6600))
# indstop = argmin(abs(kn1bit - 7400))
# cln = asarray([mean(pn1bit[indstart:indstop]), std(pn1bit[indstart:indstop])]) * 128.0
# pxrmsequiv = sqrt(cln) / (float(pixelsizearcmin * arcmins2radians))
# noisedetequiv3bit = pxrmsequiv * sqrt(float(pixelsizearcmin) / 0.0005)
#
# print(noisedetequiv3bit)
#
# plot(k, p, "g")
# #plot(kn, pn, "k")
# #plot(kn1bit, pn1bit, "r")
# #plot(kn2bit, pn2bit, "b")
# plot(kn3bit, [cln for j in pn3bit], "r")
# plot(kn3bit, pn3bit, "y")
# yscale("log")
# xscale("log")
# show()
#
#
# #-----------------------#
# #-----------------------#
# #-----------------------#







#--- READ IN INDUCED NOISE LEVELS

digitnoises1bit = zeros(len(compressions))
digitnoises1bitstd = zeros(len(compressions))
#digitnoises2bithm = zeros(len(compressions))
#digitnoises2bithmstd = zeros(len(compressions))
digitnoises2bitopt = zeros(len(compressions))
digitnoises2bitoptstd = zeros(len(compressions))
digitnoises3bit = 0
digitnoises3bitstd = 0
thmapnoise = zeros(len(compressions))

y = 0
for row in open("../fixed det noise results 1 bit/detnoiseequivalents.txt"):
    dettomapnoise = float(pixelsizearcmin) * (1.0 / (sqrt(float(pixelsizearcmin*observationlim) / raspeeds[y])))
    pxrms = noisedet / sqrt(float(pixelsizearcmin) / raspeeds[y])  # muK
    sigmamap = pxrms * float(pixelsizearcmin)  # muK arcmin
    #print(pxrms)

    rowvalues = row.split()
    # convert to map noise levels here
    digitnoises1bit[y] = float(rowvalues[0])*dettomapnoise
    digitnoises1bitstd[y] = float(rowvalues[1])*dettomapnoise
    #digitnoises2bithm[y] = float(rowvalues[2])*dettomapnoise
    #digitnoises2bithmstd[y] = float(rowvalues[3])*dettomapnoise
    digitnoises2bitopt[y] = float(rowvalues[4])*dettomapnoise
    digitnoises2bitoptstd[y] = float(rowvalues[5])*dettomapnoise
    digitnoises3bit = float(rowvalues[6]) * float(pixelsizearcmin) * (1.0 / (sqrt(float(pixelsizearcmin*128.0) / 0.0005)))
    digitnoises3bitstd = float(rowvalues[7]) * float(pixelsizearcmin) * (1.0 / (sqrt(float(pixelsizearcmin*128.0) / 0.0005)))

    pixelrms = noisedet / sqrt(float(pixelsizearcmin*observationlim) / raspeeds[y])  # muK
    thmapnoise[y] = pixelrms * float(pixelsizearcmin)  # muK arcmin
    if y < 5:
        y += 1
    else:
        pass

print(thmapnoise)
print(digitnoises1bit)


thmapnoise3bit = ( noisedet / sqrt(float(pixelsizearcmin * 128.0) / 0.0005) ) * float(pixelsizearcmin)

#--- PLOT INDUCED NOISE LEVELS
combinedcompressions = observationlim*asarray(compressions)
errorbar(combinedcompressions, digitnoises1bit-thmapnoise, yerr = digitnoises1bitstd, color = "r", marker = "o", label = "1 Bit")
#errorbar(combinedcompressions, digitnoises2bithm-thmapnoise, yerr = digitnoises2bithmstd, color = "b", marker = "o", label = "2 Bit HM")
errorbar(combinedcompressions, digitnoises2bitopt-thmapnoise, yerr = digitnoises2bitoptstd, color = "y", marker = "o", label = "2 Bit OPT")
errorbar(128.0*800000.0, digitnoises3bit-thmapnoise3bit, yerr = digitnoises3bitstd, color = "b", marker = "o", label = "3 Bit")
#plot(compressions, thmapnoise, color = "k", marker = "o", label = "Induced")
title("Induced Map Noise Level through Digitisation", fontsize = 20)
xlabel("Hits per pixel", fontsize = 20)
ylabel(r"$\Delta \sigma_{\mathrm{MAP}} \; \left[ \mu K arcmin \right]$", fontsize = 20)
xscale("log")
yscale("log")
#xlim((1, 1000000))
legend(loc = "lower left", fontsize = 15)
show()



digitnoises1bit = asarray([636.116146343, 626.724242452, 624.474736936, 619.389890122, 574.373337142, 519.37452681])
digitnoises2bithm = asarray([630.62743292, 597.00295369, 577.4905805, 558.8217488, 530.173529, 507.29169816])
digitnoises2bithmstd = asarray([170.88749743, 145.89477931, 136.36359448, 134.50887379, 125.13152214, 126.51264777])
digitnoises2bitopt = asarray([537.67434511, 532.99946591, 533.89449273, 532.98840723, 523.17172648, 512.07517825])
digitnoises2bitoptstd = asarray([152.02692924, 131.57884221, 126.07577692, 126.10295198, 123.30024347, 127.139025])

plot(compressions, digitnoises1bit-500, color = "r", marker = "o", label = "1 Bit")
errorbar(compressions, digitnoises2bithm-500, yerr = digitnoises2bithmstd, color = "b", marker = "o", label = "2 Bit HM")
errorbar(compressions, digitnoises2bitopt-500, yerr = digitnoises2bitoptstd, color = "k", marker = "o", label = "2 Bit OPT")

title("Induced Noise Level through Digitisation", fontsize = 20)
xlabel("Hits per pixel", fontsize = 20)
ylabel(r"$\sigma_{det} - 500 [\mu K \sqrt{s}]$", fontsize = 20)
xscale("log")
xlim((1, 1000000))
legend(loc = "upper left", fontsize = 15)
show()



#--- Read in CMB
cmbmap = [zeros((mappixelnumber, mappixelnumber)) for i in compressions]
cmbnoisemap2 = [zeros((mappixelnumber, mappixelnumber)) for i in compressions]
cmbnoisemap2bithm = [zeros((mappixelnumber, mappixelnumber)) for i in compressions]
cmbnoisemap2bitopt = [zeros((mappixelnumber, mappixelnumber)) for i in compressions]
cmbnoisemap1bit = [zeros((mappixelnumber, mappixelnumber)) for i in compressions]
cmbnoisemap1 = [zeros((mappixelnumber, mappixelnumber)) for i in compressions]


ci = 0
while ci < len(compressions):

    y = 0
    for row in open("../fixed det noise results 2 bit/" + str(compressions[ci]) + "hpp/cmbtemplatemapfixeddetnoise.txt"):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbmap[ci][y, x] = float(rowvalues[x])
        y += 1

    #--- Read in observed CMB
    y = 0
    for row in open("../fixed det noise results 2 bit/" + str(compressions[ci]) + "hpp/cmbnoisemapfixeddetnoise" + str(compressions[ci]) + "hppcombined.txt"):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap2[ci][y, x] += float(rowvalues[x])
        y += 1

    y = 0
    for row in open("../fixed det noise results 1 bit/" + str(compressions[ci]) + "hpp/cmbnoisemapfixeddetnoise" + str(compressions[ci]) + "hppcombined.txt"):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap1[ci][y, x] += float(rowvalues[x])
        y += 1

    #--- Read in digitised CMB
    y = 0
    for row in open("../fixed det noise results 2 bit/" + str(compressions[ci]) + "hpp/cmbnoisemap2bithmfixeddetnoise" + str(compressions[ci]) + "hppcombined.txt"):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap2bithm[ci][y, x] += float(rowvalues[x])
        y += 1

    y = 0
    for row in open("../fixed det noise results 2 bit/" + str(compressions[ci]) + "hpp/cmbnoisemap2bitoptfixeddetnoise" + str(compressions[ci]) + "hppcombined.txt"):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap2bitopt[ci][y, x] += float(rowvalues[x])
        y += 1

    y = 0
    for row in open("../fixed det noise results 1 bit/" + str(compressions[ci]) + "hpp/cmbnoisemap1bitfixeddetnoise" + str(compressions[ci]) + "hppcombined.txt"):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap1bit[ci][y, x] += float(rowvalues[x])
        y += 1

    print(ci)

    ci += 1

# --- EXTRACT POWERSPECTRUM
mapparams = [mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin]
k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap[0])
p0 = p0 * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians

detnoiseequivalents = zeros((len(compressions), 6))
mapnoiseequivalents = zeros((len(compressions), 6))

ci = 0

raspeeds = [0.0005, 0.005, 0.05, 0.5, 5.0, 50.0]
cols = [[193.0/255.0, 66.0/255.0, 66.0/255.0], [191.0/255.0, 191.0/255.0, 63.0/255.0], [63.0/255.0, 191.0/255.0, 63.0/255.0], [63.0/255.0, 63.0/255.0, 191.0/255.0], [191.0/255.0, 63.0/255.0, 127.0/255.0], [191.0/255.0, 170.0/255.0, 63.0/255.0]]




while ci < len(compressions):

    kn1, pn1, errn1, hn1 = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap1[ci])
    kn2, pn2, errn2, hn2 = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap2[ci])
    k2bithm, p2bithm, err2bithm, h2bithm = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap2bithm[ci])
    k2bitopt, p2bitopt, err2bitopt, h2bitopt = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap2bitopt[ci])
    k1bit, p1bit, err1bit, h1bit = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap1bit[ci])

    pn1 = pn1 * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    pn2 = pn2 * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    p2bithm = p2bithm * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    p2bitopt = p2bitopt * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians
    p1bit = p1bit * pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians

    # Normalise Powerspectrum
    mmi = 0
    obspow2 = 0
    obspow1 = 0
    obspow2bithm = 0
    obspow2bitopt = 0
    obspow1bit = 0

    normwindows = [[350, 1000], [350, 1000], [350, 1000], [350, 800], [50, 300], [15, 100]]

    while mmi < len(k0):
        if k0[mmi] >= normwindows[ci][0] and k0[mmi] <= normwindows[ci][1]:
            obspow2 += pn2[mmi] ** 2.0
            obspow1 += pn1[mmi] ** 2.0
            obspow2bithm += p2bithm[mmi] ** 2.0
            obspow2bitopt += p2bitopt[mmi] ** 2.0
            obspow1bit += p1bit[mmi] ** 2.0

        mmi += 1
    p2bithm = p2bithm * (obspow2 / obspow2bithm) ** 0.5
    p2bitopt = p2bitopt * (obspow2 / obspow2bitopt) ** 0.5
    p1bit = p1bit * (obspow1 / obspow1bit) ** 0.5


    plot(kn2, pn2, color = tuple(cols[ci] + [1.0]), label = str(compressions[ci])+" hpp")
    #plot(k2bithm, p2bithm, color = tuple(cols[ci] + [0.7]), ls = "--")
    #plot(k2bitopt, p2bitopt, color = tuple(cols[ci] + [0.7]), ls = "-.")
    plot(k1bit, p1bit, color = tuple(cols[ci] + [0.7]), ls = ":")


    # --- DETERMINE INDUCED NOISE LEVEL
    k1bitflatstart = [5900, 4500, 3500, 2500, 2000, 1500]
    k1bitflatstop = 6780
    k2bitflatstart = [5900, 4500, 3500, 3000, 2000, 1500]
    k2bitflatstop = 6780
    ind1bitstart = argmin(abs(k1bit - k1bitflatstart[ci]))
    ind1bitstop = argmin(abs(k1bit - k1bitflatstop))
    ind2bitstart = argmin(abs(k2bithm - k2bitflatstart[ci]))
    ind2bitstop = argmin(abs(k2bithm - k2bitflatstop))
    noisecl2bithmcombined = asarray([mean(p2bithm[ind2bitstart:ind2bitstop]), std(p2bithm[ind2bitstart:ind2bitstop])])
    noisecl2bitoptcombined = asarray([mean(p2bitopt[ind2bitstart:ind2bitstop]), std(p2bitopt[ind2bitstart:ind2bitstop])])
    noisecl1bitcombined = asarray([mean(p1bit[ind1bitstart:ind1bitstop]), std(p1bit[ind1bitstart:ind1bitstop])])


    #print(noisecl2bithmcombined)

    noisecl2bithm = noisecl2bithmcombined * float(observationlim)  # muK^2
    noisecl2bitopt = noisecl2bitoptcombined * float(observationlim)  # muK^2
    noisecl1bit = noisecl1bitcombined * float(observationlim)  # muK^2
    pixelrmsequiv2bithm = sqrt(noisecl2bithm) / (float(pixelsizearcmin * arcmins2radians))  # muK
    pixelrmsequiv2bitopt = sqrt(noisecl2bitopt) / (float(pixelsizearcmin * arcmins2radians))  # muK
    pixelrmsequiv1bit = sqrt(noisecl1bit) / (float(pixelsizearcmin * arcmins2radians))  # muK
    noisedetequiv2bithm = pixelrmsequiv2bithm * sqrt(float(pixelsizearcmin) / raspeeds[ci])  # muK sqrt(s)
    noisedetequiv2bitopt = pixelrmsequiv2bitopt * sqrt(float(pixelsizearcmin) / raspeeds[ci])  # muK sqrt(s)
    noisedetequiv1bit = pixelrmsequiv1bit * sqrt(float(pixelsizearcmin) / raspeeds[ci])  # muK sqrt(s)
    # noisemapequiv1bit = pixelrmsequiv1bit * float(pixelsizearcmin)  # muK arcmin
    # noisemapequiv2bithm = pixelrmsequiv2bithm * float(pixelsizearcmin)  # muK arcmin
    # noisemapequiv2bitopt = pixelrmsequiv1bit * float(pixelsizearcmin)  # muK arcmin


    for j in range(2):
        detnoiseequivalents[ci][0+j] = noisedetequiv1bit[j]
        detnoiseequivalents[ci][2+j] = noisedetequiv2bithm[j]
        detnoiseequivalents[ci][4+j] = noisedetequiv2bitopt[j]

    print("--------------------")
    print("HPP " + str(compressions[ci]))
    print("DETECTOR NOISE EQUIVALENT")
    print("2 BIT HM:")
    print(noisedetequiv2bithm)
    print("2BIT OPT:")
    print(noisedetequiv2bitopt)
    print("--------------------")

    print(ci)

    ci += 1

#--- SAVE NOISE LVL EQUIVALENTS
savetxt("detnoiseequivalents.txt", detnoiseequivalents)


#--- PLOT POWERSPECTRUM
plot(k0, p0, "k", label = "True")
legend(loc="lower left", fontsize = 15)
title(r"CMB Powerspectrum 2 Bit Digitisation", fontsize = 20)
ylabel(r"$C_l [\mu K^2]$", fontsize = 20)
xlabel(r"Multipole Moment $l$", fontsize = 20)
xscale("log")
yscale("log")
show()

#kresid, presid, errresid, hresid = cosm.sriniPowerSpectrum(mapparams, noiseresid)
#kdresid, pdresid, errdresid, hdresid = cosm.sriniPowerSpectrum(mapparams, digitnoiseresid)
#k2bithm, p2bithm, err2bithm, h2bithm = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap2bithm)
#k2bitopt, p2bitopt, err2bitopt, h2bitopt = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap2bitopt)


# Normalise Powerspectrum
mmi = 0
#truepow = 0
obspow1 = 0
#obspow2 = 0
obspow1bit = 0
#obspow2bithm = 0
#obspow2bitopt = 0
while mmi < len(k0):
    if k0[mmi] >= 350 and k0[mmi] <= 1000:
        #truepow += p0[mmi]**2.0
        obspow1 += pn[mmi]**2.0
        #obspow2 += pn2[mmi]**2.0
        obspow1bit += p1bit[mmi]**2.0
        #obspow2bithm += p2bithm[mmi]**2.0
        #obspow2bitopt += p2bitopt[mmi]**2.0

    mmi += 1

#p = p * (truepow/obspow)**0.5
p1bit = p1bit * (obspow1/obspow1bit)**0.5
#p2bithm = p2bithm * (obspow2/obspow2bithm)**0.5
#p2bitopt = p2bitopt * (obspow2/obspow2bitopt)**0.5

#--- DETERMINE INDUCED NOISE LEVEL
k1bitflat = 5900
ind = argmin(abs(k1bit-k1bitflat))
noisecl1bitcombined = mean(p1bit[ind:])
noisecl1bit = noisecl1bitcombined*float(observationlim) # muK^2
pixelrmsequiv = sqrt(noisecl1bit)/(float(pixelsizearcmin*arcmins2radians)) # muK
noisedetequiv = pixelrmsequiv*sqrt(float(pixelsizearcmin)/raspeed) # muK sqrt(s)
noisemapequiv = pixelrmsequiv*float(pixelsizearcmin) # muK arcmin

print("DETECTOR NOISE EQUIVALENT")
print(noisedetequiv)


#--- CALCULATE DL PS

dl0 = p0 * k0 * (k0 + 1.0) / (2.0 * pi)
dln = pn * kn * (kn + 1.0) / (2.0 * pi)
dl1bit = p1bit * k1bit * (k1bit + 1.0) / (2.0 * pi)
#dlobstot2bithm = float(pixelnumber) * float(pixelnumber) * p2bithm * k2bithm * (k2bithm + 1.0) / pi
#dlobstot2bitopt = float(pixelnumber) * float(pixelnumber) * p2bitopt * k2bitopt * (k2bitopt + 1.0) / pi


#--- PLOT PS

subplot(1, 2, 1)
plot(k0, p0, "b", label="True")
plot(kn, pn, "g", label="Observed")
plot(k1bit, p1bit, "r", label="Digitised")
plot(k0, [noiseclcombined for l in k0], "k", label="Noise Induced")
plot(k0, [noisecl1bitcombined for u in k0], "k", ls = "dashed", label = "Noise Digitisation")
title(r"CMB Powerspectrum Fixed Detector Noise", fontsize = 20)
xlabel(r"Mutlipole Moment $l$", fontsize = 20)
ylabel(r"$C_l [\mu K^2]$", fontsize = 20)
xscale("log")
yscale("log")
legend(loc = "lower left", fontsize = 15)


subplot(1, 2, 2)
plot(k0, p0, "b", label="True")
plot(kn, pn, "g", label="Observed")
plot(k1bit, p1bit, "r", label="Digitised")
plot(k0, [noiseclcombined for l in k0], "k", label="Noise")
plot(k0, [noisecl1bitcombined for u in k0], "k", ls = "dashed")
title(r"Zoomed In", fontsize = 20)
xlabel(r"Mutlipole Moment $l$", fontsize = 20)
ylabel(r"$C_l [\mu K^2]$", fontsize = 20)
xscale("log")
yscale("log")
xlim((3000, max(k0)))
ylim((1e-8, 1e-4))

show()