# Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from scipy.stats import binned_statistic
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from cosmdigitclasses import *
from numpy import *
from scipy.signal import convolve2d
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
    kx = 2*pi*fftshift(fftfreq(pixelnumber, pixelspacingrad))
    kx2d = tile(kx, (pixelnumber, pixelnumber))
    ky = 2*pi*fftshift(fftfreq(pixelnumber, pixelspacingrad))
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
cmbfreqspace = (realpart + 1j*imagpart)

# Transform into map
cmbmap = fft.ifft2(cmbfreqspace)[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
cmbmap = real(cmbmap)

k, p_k, p_kerr, hitcount = cosm.sriniPowerSpectrum([512, 512, 2, 2], cmbmap)
k = asarray(k)
p_k = asarray(p_k)
c_l = p_k*k*k + p_k*k
hitcount = asarray(hitcount)
plot(k, c_l)
title("Powerspectrum")
xlabel("Multipole moment k")
ylabel("p_k*k(k+1)")
show()

dasd321



cmbnoisemap = convolve2d(cmbmap, normal(size = shape(cmbmap)))[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
imshow(cmbnoisemap)
title("CMB and Noise by convolution")
colorbar()
show()


k, p_k, p_kerr, hitcount = cosm.sriniPowerSpectrum([512, 512, 2, 2], cmbmap)
k = asarray(k)
hitcount = asarray(hitcount)
plot(k, p_k)
show()

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
        # Create noise
        noisereal = normal(size = compression, scale = 0.0006)
        noiseimag = normal(size = compression, scale = 0.0006)
        noisefreqspace = noisereal + 1j*noiseimag
        noise = real(ifft(noisefreqspace))
        # Add noise and signal as time stream
        tod = cmbmap[d, ri] + noise
        cesscans[d, rstart:rstop] = tod

print("Added noise")

# Recompress into map
cmbnoisemap = zeros((nodecscans, norablocks))
for d in range(shape(cesscans)[0]):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        cmbnoisemap[d, ri] = mean(cesscans[d, rstart:rstop])

print("Recompressed")
imshow(cmbnoisemap)
show()

mapparams = [512, 512, 2, 2]
k, p, err, h = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap)
k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
plot(k, p, "r")
plot(k0, p0, "b")
diff = asarray(k)-asarray(k0)
xlabel("l")
ylabel("cl")
xscale("log")
yscale("log")
show()

plot(k, diff)
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