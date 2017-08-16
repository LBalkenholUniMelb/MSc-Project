# Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from scipy.stats import binned_statistic
from numpy.fft import ifft2, fft2
from cosmdigitclasses import *
from numpy import *
rc('text', usetex=True)

cosm = Cosmologist()

# Define map parameters
fieldsizearcmins = 2048
pixelsizearcmin = 2
pixelnumber = 1024
df = 1

# Read in data
filelocation = "for_lennart/plik_plus_r0p01_highell_lensedtotCls.dat"
file = open(filelocation)
l = []
dl = []
for row in file:
    rowvalues = [float(i) for i in row.split()]
    l.append(rowvalues[0])
    dl.append(rowvalues[1])

# Calculate cl
cl = [(dl[i]*2*pi)/(l[i]*(l[i]+1)) for i in range(len(l))]

# Create 2dCl
cl2d = zeros((pixelnumber, pixelnumber))
readincld2 = True
if readincld2:
    y = 0
    for row in open("cl2d.txt"):
        rowvalues = row.split()
        for x in range(pixelnumber):
            cl2d[y, x] = float(rowvalues[x])
        y += 1
else:
    pixelspacingrad = (pixelsizearcmin/60)*(pi/180)
    kx = 2*pi*fftfreq(pixelnumber, pixelspacingrad)
    kx2d = tile(kx, (pixelnumber, pixelnumber))
    ky = 2*pi*fftfreq(pixelnumber, pixelspacingrad)
    ky2d = transpose(tile(ky, (pixelnumber, pixelnumber)))
    cl2d = zeros((pixelnumber, pixelnumber))

    for y in range(pixelnumber):
        for x in range(pixelnumber):
            k = sqrt(kx2d[y, x]*kx2d[y, x] + ky2d[y, x]*ky2d[y, x])
            ind = argmin([abs(i-k) for i in l])
            cl2d[y, x] = cl[ind]

# Combine noise and cmb element-wise
factor = sqrt((df/2)*cl2d)

realpart = factor * normal(size = (pixelnumber, pixelnumber))
imagpart = factor * normal(size = (pixelnumber, pixelnumber))
cmbnoisefreqspace = (realpart + 1j*imagpart)

# Transform into map
cmbnoisemap = fft.ifft2(cmbnoisefreqspace)[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
cmbnoisemap = real(cmbnoisemap)
psvmap = sum(cmbnoisemap*cmbnoisemap)

# k, p_k, p_kerr, hitcount = cosm.sriniPowerSpectrum([512, 512, 2, 2], cmbnoisemap)
# p_k = asarray(p_k)
# hitcount = asarray(hitcount)
# psvrad = sum(p_k*hitcount)

# m = mean(dpsv)
# s = std(dpsv)
# plot(psvindex, dpsv)
# xscale("linear")
# yscale("linear")
# title(r"$\Delta P$ for 100 Realisations")
# xlabel(r"Iteration number")
# ylabel(r"$\Delta P$")
# plot(psvindex, [m for i in psvindex], "0.0")
# plot(psvindex, [m+s for i in psvindex], "0.5")
# plot(psvindex, [m-s for i in psvindex], "0.5")
# show()

# Recreate Scan Strategy
declims = [0, 1024] #arcmins
ralims = [0, 1024] #arcmins
readoutfreq = 6 #Hz
raspeed = 0.1 #arcmin/s
nodecscans = 512
norablocks = 512
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
cesscans = zeros((nodecscans, radatapoints))

for d in range(nodecscans):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        # Create relevant tod signals here
        # sig = cmbnoisemap[d, ri]*normal(2.7251, 0.0006, compression)
        cesscans[d, rstart:rstop] = asarray([cmbnoisemap[d, ri]] * compression)

# Recompress into map
maprecreated = zeros((nodecscans, norablocks))
for d in range(shape(cesscans)[0]):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        m = mean(cesscans[d, rstart:rstop])
        maprecreated[d, ri] = m

# Get power spectrum
k, p_k, p_kerr, hitcount = cosm.sriniPowerSpectrum([512, 512, 2, 2], maprecreated)
k0, p_k0, p_kerr0, hitcount0 = cosm.sriniPowerSpectrum([512, 512, 2, 2], cmbnoisemap)
p_k = asarray(p_k)
k = asarray(k)
c_l = (k*k*p_k + k*p_k)/(2*pi)
p_k0 = asarray(p_k0)
k0 = asarray(k0)
c_l0 = (k0*k0*p_k0 + k0*p_k0)/(2*pi)
hitcount = asarray(hitcount)
plot(k, c_l-c_l0, "r")
#plot(k0, c_l0, "b")
xscale("log")
yscale("linear")
show()

# # Compare powerspectra
# mapparams = [512, 512, 2, 2]
# kfull, clfull = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap, col = "b")
# kcomp, clcomp = cosm.sriniPowerSpectrum(mapparams, maprecreated, col = "None", mark="o")
# plot(kfull, clfull, "b", label="Orgininal")
# plot(kcomp, clcomp, "r", label="Recreated")
# xlabel("l")
# ylabel("cl")
# xscale("log")
# yscale("log")
# title("Powerspectrum of CMB and Noise Map, original and recreated")
# legend()
# show()