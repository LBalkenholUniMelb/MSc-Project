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

var = 6*10**(-6)#std(cmbmap)
avg = 0

# Recreate Scan Strategy
declims = [0, 1024] #arcmins
ralims = [0, 1024] #arcmins
readoutfreq = 6 #Hz
raspeed = 0.1 #arcmin/s
nodecscans = 512
norablocks = 512
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)

observationno = range(1)
observations = [zeros((nodecscans, norablocks)) for i in observationno]
cesscans = [zeros((nodecscans, radatapoints)) for i in observationno]

# Decompose into bolometer signals

for d in range(nodecscans):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        # Create and noise
        noise = normal(avg, var, size = (len(observationno), compression))
        for obs in observationno:
            tod = cmbmap[d, ri] + noise[obs] # digitise here
            tod = digitise1bit(tod, var)
            cesscans[obs][d, rstart:rstop] = tod

    print(d)

# Recompress into map

print("Decomposition completed")


for d in range(shape(cesscans[0])[0]):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        for obs in observationno:
            observations[obs][d, ri] = mean(cesscans[obs][d, rstart:rstop])

    print(d)

print("Beginning PS extraction")

cmbnoisemap = zeros((nodecscans, norablocks))
for obs in observations:
    cmbnoisemap = cmbnoisemap + obs
cmbnoisemap = cmbnoisemap * (1.0/len(observationno))

# Plot powerspectrum
mapparams = [512, 512, 2, 2]
k, p, err, h = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap)
k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
subplot(1, 2, 1)
plot(k, p, "r", label = "CMB+N")
plot(k0, p0, "b", label = "CMB")
title(r"Powerspectrum")
xlabel(r"Wavevector $k$")
ylabel(r"Power $P_k$")
xscale("log")
yscale("log")
legend()
subplot(1, 2, 2)
plot(k0, asarray(p0)-asarray(p))
title(r"Powerspectrum difference")
xlabel(r"Wavevector $k$")
ylabel(r"Difference in $P_k$")
xscale("log")
print(mean(asarray(p0)-asarray(p)))
print(std(asarray(p0)-asarray(p)))
show()