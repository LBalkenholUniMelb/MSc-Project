#--- Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from cosmdigitclasses import *
from numpy import *
from scipy.signal import convolve2d
from scipy.stats import *


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
readoutfreq = 400 #Hz
raspeed = 0.001 #arcmin/s
nodecscans = mappixelnumber
norablocks = mappixelnumber
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
observationlim = 1

#--- Define noise level
noisemap = 3.0 # muK arcmin
noisepix = noisemap/float(pixelsizearcmin)
fsky = (mapfieldsize*mapfieldsize)/(4.0*pi*60.0*60.0*(180.0/pi)**2.0)
noisecl = 4.0*pi*fsky*noisepix*noisepix/float(mappixelnumber*mappixelnumber)
eta = sqrt(float(compression)) * sqrt(float(observationlim)) * sqrt(float(noisecl) / (float(pixelnumber * pixelnumber)))

#--- Read in true CMB
cmbmap = zeros((mappixelnumber, mappixelnumber))
y = 0
for row in open("results/cmbmaparblvltemplate.txt"):
    rowvalues = row.split()
    for x in range(mappixelnumber):
        cmbmap[y, x] = float(rowvalues[x])
    y += 1

#--- Read in CMB Observations
cmbnoisemap2 = zeros((mappixelnumber, mappixelnumber))
cmbnoisemap1 = zeros((mappixelnumber, mappixelnumber))
cmbnoisemap1bit = zeros((mappixelnumber, mappixelnumber))
cmbnoisemap2bithm = zeros((mappixelnumber, mappixelnumber))
cmbnoisemap2bitopt = zeros((mappixelnumber, mappixelnumber))

ranks = 10
for rank in range(ranks):

    noisemap2name = "results/cmbnoisemap" + str(rank+10) + ".txt"
    y = 0
    for row in open(noisemap2name):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap2[y, x] += float(rowvalues[x])
        y += 1

    noisemap1name = "results/cmbnoisemap" + str(rank) + ".txt"
    y = 0
    for row in open(noisemap1name):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap1[y, x] += float(rowvalues[x])
        y += 1

    digit1bitnoisemapname = "results/cmbnoisemap1bit" + str(rank) + ".txt"
    y = 0
    for row in open(digit1bitnoisemapname):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap1bit[y, x] += float(rowvalues[x])
        y += 1

    digithmnoisemapname = "results/cmbnoisemap2bithm" + str(rank+10) + ".txt"
    y = 0
    for row in open(digithmnoisemapname):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap2bithm[y, x] += float(rowvalues[x])
        y += 1

    digitoptnoisemapname = "results/cmbnoisemap2bitopt" + str(rank+10) + ".txt"
    y = 0
    for row in open(digitoptnoisemapname):
        rowvalues = row.split()
        for x in range(mappixelnumber):
            cmbnoisemap2bitopt[y, x] += float(rowvalues[x])
        y += 1

    print("COMPLETED READING IN " + str(100*rank/ranks) + "%")

cmbnoisemap1 = cmbnoisemap1 / float(ranks)
cmbnoisemap2 = cmbnoisemap2 / float(ranks)
cmbnoisemap1bit = cmbnoisemap1bit / float(ranks)
cmbnoisemap2bithm = cmbnoisemap2bithm / float(ranks)
cmbnoisemap2bitopt = cmbnoisemap2bitopt / float(ranks)

#noiseresid = cmbnoisemap-cmbmap
#digitnoiseresid = cmbnoisemap1bit-sqrt(sum(cmbnoisemap1bit**2.0)/sum(cmbnoisemap**2.0))*cmbnoisemap

#--- THEORY PLOTS INTERRUPT


cmbmap1adj = cmbmap*sqrt(sum(cmbnoisemap1**2.0)/sum(cmbmap**2.0))
cmbnoisemap1bitadj = cmbnoisemap1bit*sqrt(sum(cmbnoisemap1**2.0)/sum(cmbnoisemap1bit**2.0))

cmbmap2adj = cmbmap*sqrt(sum(cmbnoisemap2**2.0)/sum(cmbmap**2.0))
cmbnoisemap2bithmadj = cmbnoisemap2bithm*sqrt(sum(cmbnoisemap2**2.0)/sum(cmbnoisemap2bithm**2.0))
cmbnoisemap2bitoptadj = cmbnoisemap2bitopt*sqrt(sum(cmbnoisemap2**2.0)/sum(cmbnoisemap2bitopt**2.0))


subplot(1, 2, 1)
imshow(cmbnoisemap1bitadj)
title("CMB Map Digitised Observation", fontsize = 20)
colorbar()

subplot(1, 2, 2)
imshow(cmbnoisemap1bitadj-cmbnoisemap1)
title("CMB Map Digitised - CMB Map Observed", fontsize = 20)
colorbar()

show()

plot(cmbmap2adj.flatten(), (cmbmap2adj-cmbnoisemap2bithmadj).flatten(), "bx", label = "2 Bit HM")
plot(cmbmap2adj.flatten(), (cmbmap2adj-cmbnoisemap2bitoptadj).flatten(), "rx", label = "2 Bit OPT")
plot(cmbmap1adj.flatten(), (cmbmap1adj-cmbnoisemap1bitadj).flatten(), "kx", label = "1 Bit")
theo = linspace(-4, 4, 100)
plot(theo, theo-(norm.cdf(theo)-norm.cdf(-theo)), "y", label = "Theoretical")
legend()
xlim((-0.0002, 0.0002))
ylim((-0.0002, 0.0002))
show()

#--- CONTINUE

mapparams = [mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin]
k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
kn1, pn1, errn1, hn1 = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap1)
kn2, pn2, errn2, hn2 = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap2)
k1bit, p1bit, err1bit, h1bit = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap1bit)
#kresid, presid, errresid, hresid = cosm.sriniPowerSpectrum(mapparams, noiseresid)
#kdresid, pdresid, errdresid, hdresid = cosm.sriniPowerSpectrum(mapparams, digitnoiseresid)
k2bithm, p2bithm, err2bithm, h2bithm = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap2bithm)
k2bitopt, p2bitopt, err2bitopt, h2bitopt = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap2bitopt)


# Normalise Powerspectrum
mmi = 0
truepow = 0
obspow1 = 0
obspow2 = 0
obspow1bit = 0
obspow2bithm = 0
obspow2bitopt = 0
while mmi < len(k0):
    if k0[mmi] >= 350 and k0[mmi] <= 1000:
        truepow += p0[mmi]**2.0
        obspow1 += pn1[mmi]**2.0
        obspow2 += pn2[mmi]**2.0
        obspow1bit += p1bit[mmi]**2.0
        obspow2bithm += p2bithm[mmi]**2.0
        obspow2bitopt += p2bitopt[mmi]**2.0
    mmi += 1

#p = p * (truepow/obspow)**0.5
p1bit = p1bit * (obspow1/obspow1bit)**0.5
p2bithm = p2bithm * (obspow2/obspow2bithm)**0.5
p2bitopt = p2bitopt * (obspow2/obspow2bitopt)**0.5

p = pn1
k = kn1

#--- PLOT RESIDUAL VAL PLOTS

hmresid = cmbnoisemap2bithm-sqrt(sum(cmbnoisemap2bithm**2.0)/sum(cmbnoisemap2**2.0))*cmbnoisemap2
optresid = cmbnoisemap2bitopt-sqrt(sum(cmbnoisemap2bitopt**2.0)/sum(cmbnoisemap2**2.0))*cmbnoisemap2

# digitnoise1 = cmbnoisemap1bit-sqrt(sum(cmbnoisemap1bit**2.0)/sum(cmbnoisemap1**2.0))*cmbnoisemap1
#
# subplot(1, 2, 1)
# imshow(cmbnoisemap1*sqrt((sum(cmbnoisemap1bit**2.0))/(sum(cmbnoisemap1**2.0))))
# colorbar()
# title("CMB", fontsize = 20)
#
# subplot(1, 2, 2)
# imshow(cmbnoisemap1bit-sqrt(sum(cmbnoisemap1bit**2.0)/sum(cmbnoisemap1**2.0))*cmbnoisemap1)
# colorbar()
# title("CMB-POW*1BIT", fontsize = 20)
#
# show()


#--- PLOT MAPS
#
# subplot(2, 2, 1)
# imshow(cmbmap)
# #colorbar()
# title("CMB", fontsize = 20)
#
# subplot(2, 2, 2)
# imshow(cmbnoisemap)
# #colorbar()
# title("CMB+N", fontsize = 20)
#
# subplot(2, 2, 3)
# imshow(cmbnoisemap2bithm)
# #colorbar()
# title("2 BIT HM", fontsize = 20)
#
# subplot(2, 2, 4)
# imshow(cmbnoisemap2bitopt)
# #colorbar()
# title("2 BIT OPT", fontsize = 20)
#
#
# show()
#
#
#
#--- PLOT DL PS

dltot = float(pixelnumber) * float(pixelnumber) * p0 * k0 * (k0 + 1.0) / pi
dlobstot = float(pixelnumber) * float(pixelnumber) * p * k * (k + 1.0) / pi
dlobstot1bit = float(pixelnumber) * float(pixelnumber) * p1bit * k1bit * (k1bit + 1.0) / pi
dlobstot2bithm = float(pixelnumber) * float(pixelnumber) * p2bithm * k2bithm * (k2bithm + 1.0) / pi
dlobstot2bitopt = float(pixelnumber) * float(pixelnumber) * p2bitopt * k2bitopt * (k2bitopt + 1.0) / pi

cltot = (dltot*pi)/(k0*(k0+1.0))
clobstot = (dlobstot*pi)/(k*(k+1.0))
#clobstot1bit = (dlobstot1bit*pi)/(k1bit*(k1bit+1.0))
clobstot2bithm = (dlobstot2bithm*pi)/(k2bithm*(k2bithm+1.0))
clobstot2bitopt = (dlobstot2bitopt*pi)/(k2bitopt*(k2bitopt+1.0))

subplot(1, 2, 1)
plot(k0, dltot, label="True")
plot(k, dlobstot, label="Observed")
plot(k1bit, dlobstot1bit, label="1 Bit")
#plot(k2bitopt, dlobstot2bitopt, label="OPT")
legend(fontsize = 15)
xscale("linear")
yscale("linear")
ylabel(r"$D_l [\mu K^2]$", fontsize = 20)
xlabel(r"Multipole Moment $l$", fontsize = 20)
xlim((0, 2500))
title(r"Powerspectrum", fontsize = 20)

subplot(1, 2, 2)
plot(k1bit, dlobstot1bit-dlobstot)
#plot(k0, dlobstot2bitopt-dlobstot, label="OPT")
title(r"Digitised - Observed PS", fontsize = 20)
xlabel(r"Multipole Moment $l$", fontsize = 20)
ylabel(r"$D_l^{\mathrm{Digit}} - D_l^{\mathrm{Obs}}$", fontsize = 20)
xscale("linear")
yscale("linear")
#ylim((0, 9))
legend()

show()


#--- PLOT CL PS

subplot(1, 2, 1)
plot(k0, cltot, label="True")
plot(k, clobstot, label="Observed")
plot(k2bithm, clobstot2bithm, label="HM")
#plot(k2bitopt, clobstot2bitopt, label="OPT")
legend(fontsize = 15, loc = 3)
xscale("log")
yscale("log")
ylabel(r"$C_l [\mu K^2]$", fontsize = 20)
xlabel(r"Multipole Moment $l$", fontsize = 20)
#xlim((0, 2500))
title(r"Powerspectrum", fontsize = 20)

subplot(1, 2, 2)
plot(k0, cltot, label="True")
plot(k, clobstot, label="Observed")
plot(k2bithm, clobstot2bithm, label="HM")
#plot(k2bitopt, clobstot2bitopt, label="OPT")
#legend(fontsize = 15)
xscale("log")
yscale("log")
#ylabel(r"$C_l [\mu K^2]$", fontsize = 20)
xlabel(r"Multipole Moment $l$", fontsize = 20)
xlim((3e3, 8e3))
ylim((1e-8, 1e-5))
title(r"Zoomed in", fontsize = 20)

show()


# #--- PLOT DRESID PS
# plot(kdresid, pdresid)
# title("CMBND-CMBN PS")
# xscale("log")
# yscale("log")
# xlim((min(kdresid), max(kdresid)))
# show()