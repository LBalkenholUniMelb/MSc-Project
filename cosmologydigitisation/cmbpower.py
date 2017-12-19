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
df = 1.0/(float(fieldsizearcmins)*arcmins2radians)

# Read in data
filelocation = "for_lennart/plik_plus_r0p01_highell_lensedtotCls.dat"
file = open(filelocation)
l = []
dl = []
for row in file:
    rowvalues = [float(i) for i in row.split()]
    l.append(rowvalues[0])
    dl.append(rowvalues[1])

dlinput = asarray(dl)
linput = asarray(l)

# Calculate cl
cl = [(dl[i]*2.0*pi)/(l[i]*(l[i]+1.0)) for i in range(len(l))]

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
    pixelspacingrad = (float(pixelsizearcmin)/60.0)*(pi/180.0)
    kx = 2.0*pi*fftshift(fftfreq(pixelnumber, pixelspacingrad))
    kx2d = tile(kx, (pixelnumber, pixelnumber))
    ky = 2.0*pi*fftshift(fftfreq(pixelnumber, pixelspacingrad))
    ky2d = transpose(tile(ky, (pixelnumber, pixelnumber)))
    cl2d = zeros((pixelnumber, pixelnumber))

    for y in range(pixelnumber):
        for x in range(pixelnumber):
            k = sqrt(kx2d[y, x]*kx2d[y, x] + ky2d[y, x]*ky2d[y, x])
            ind = argmin([abs(i-k) for i in l])
            cl2d[y, x] = cl[ind]
        print(y)
    savetxt("cl2drevised.txt", cl2d)

#print("POWER IN DL:")
#print(sum(dl*conjugate(dl)))

#print("POWER IN CL2D:")
#print(sum(cl2d*conjugate(cl2d)))


# Combine noise and cmb element-wise
factor = sqrt((df/2.0)*cl2d)

print("POW IN FAC:")
print(sum(factor*factor))

#print("POWER IN FACTOR")
#print(sum(factor*conjugate(factor)))

# file = open("cmbnoiseseed.txt")
# re = zeros((pixelnumber, pixelnumber))
# im = zeros((pixelnumber, pixelnumber))
# rowindex = 0
# for row in file:
#     rowvalues = [float(i) for i in row.split()]
#     re[rowindex] = rowvalues[:pixelnumber]
#     im[rowindex] = rowvalues[pixelnumber:]
#     rowindex += 1

cmbmapavg = zeros((pixelnumber/2, pixelnumber/2))

i = 0
dltot = 0
first = True
powdiffs = []
powfreqs = []
noreals = 500.0
while i < int(noreals):

    re = normal(0, 1, (pixelnumber, pixelnumber))
    im = normal(0, 1, (pixelnumber, pixelnumber))

    realpart = factor * re
    imagpart = factor * im
    cmbfreqspace = (realpart + 1.0j*imagpart)

    powfreq = sum(cmbfreqspace*conjugate(cmbfreqspace))
    powfreqs.append(powfreq)

    # Transform into map
    cmbmap = fft.ifft2(cmbfreqspace)[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
    #print("POWER IN CMB REAL SPACE * 4 * N^2:")
    #print(sum(cmbmap*conjugate(cmbmap))*4.0*pixelnumber*pixelnumber)
    #print(sum(cmbmap*conjugate(cmbmap)))

    cmbmap = real(cmbmap)

    cmbmapavg += cmbmap*4.0*pixelnumber*pixelnumber


    #print("BECAUSE WE DISCARD IMAGINARY PART:")
    #print(sum(cmbmap*cmbmap*2.0))
    #print(sum(cmbmap*cmbmap))

    k, p, perr, h = cosm.sriniPowerSpectrum([512, 512, 2, 2], cmbmap)

    #print("ULTIMATE:")
    powps = sum(p*h)*8.0*1024.0*1024.0
    powfac = (4.0*1024.0*1024.0)


    powdiffs.append((powfreq-powps)/powfreq)

    dl = p*k*(k+1.0)/(2.0*pi)
    dl = dl*powfac

    #dl = dl*4.0e6

    if first:
        dltot = dl
        first = False
    else:
        dltot = dltot + dl
    print(i)
    i += 1

means = []
j = 0
while j < i:
    l = 0
    m = 0
    while l <= j:
        m += powdiffs[l]
        l += 1
    means.append(m/(j+1.0))
    j += 1

xscale("linear")
yscale("linear")
cmbmapavg = cmbmapavg/noreals
imshow(cmbmapavg)
colorbar()
show()

xscale("linear")
yscale("linear")
powfreqs = asarray(powfreqs)/2.0
m = mean(powfreqs)
s = std(powfreqs)
plot(range(i), powfreqs, "k")
plot(range(i), [m for l in range(i)], "b")
plot(range(i), [m+s for l in range(i)], "b")
plot(range(i), [m-s for l in range(i)], "b")
plot(range(i), [sum(factor*factor) for p in range(i)], "r")
xscale("linear")
yscale("linear")
show()


plot(range(i), means)
xscale("linear")
yscale("linear")
show()

plot(range(i), powdiffs)
xscale("linear")
yscale("linear")
show()

subplot(2, 1, 1)
title("Powerspectrum", fontsize = 25)
dltot = dltot/noreals
plot(k, dltot, "r", label = "Reconstructed")
plot(linput, dlinput, "b", label= "Input")
xlabel(r"Multipole Moment $l$", fontsize = 20)
ylabel(r"$l(l+1)C_l/2\pi [\mu \mathrm{K}^2]$", fontsize = 20)
legend(fontsize = 20)
xscale("linear")
yscale("linear")
xlim((0, 2000))

# find where linput and k coincide
keq = []
psdiff = []
for count, wavevec in enumerate(k):
    ind = argmin(abs(linput-wavevec))
    keq.append(linput[ind])
    psdiff.append(abs(dlinput[ind]-dltot[count]))

subplot(2, 1, 2)
plot(keq, psdiff)
title("Powerspectrum difference", fontsize = 25)
xlabel(r"Multipole Moment $l$", fontsize = 20)
ylabel(r"Difference in $[\mu \mathrm{K}^2]$", fontsize = 20)
xscale("linear")
yscale("linear")
xlim((0, 2000))

show()



#print("POWER IN CMB REAL SPACE:")
#print(sum(cmbmap*conjugate(cmbmap))*2)
#print(sum(cmbmap*conjugate(cmbmap)))

#print("POWER IN P:")
#print(hits)
#print(sum(p*hitsc))

#plot(l, dl, "b")

# Recreate Scan Strategy
declims = [0, 1024] #arcmins
ralims = [0, 1024] #arcmins
readoutfreq = 6 #Hz
raspeed = 0.1 #arcmin/s
nodecscans = 512
norablocks = 512
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)

noisemap = 2.5 #muK
noisepixel = noisemap/sqrt(fieldsizearcmins*fieldsizearcmins*0.5*0.5) #muK/arcmin
noisetod = noisepixel*sqrt(float(compression))

observationno = range(1)
observations = [zeros((nodecscans, norablocks)) for i in observationno]
cesscans = [zeros((nodecscans, radatapoints)) for i in observationno]

# Decompose into bolometer signals
noise = normal(0, noisetod, size = (len(observationno), nodecscans, norablocks*compression))

for d in range(nodecscans):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        # Create and noise
        for obs in observationno:
            tod = cmbmap[d, ri] + noise[obs][d][rstart:rstop]
            #print(std(tod))
            # adding noise to cmb map
            #cmbmap[d, ri] = mean(tod)
            # digitise here
            todpow = sum(asarray(tod)*asarray(tod))
            digitise1bit(tod, 1)
            #toddigitpow = sum(asarray(tod)*asarray(tod))
            #print(((todpow/toddigitpow)**0.5))
            #tod = ((todpow/toddigitpow)**0.5)*tod
            cesscans[obs][d, rstart:rstop] = tod

    print("Noise & Digitisation: " + str(int(100*d/nodecscans)) + "%")

# Recompress into map

print("Decomposition completed")


for d in range(shape(cesscans[0])[0]):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        for obs in observationno:
            observations[obs][d, ri] = mean(cesscans[obs][d, rstart:rstop])

    print("Compression: " + str(int(100*d/nodecscans)) + "%")

print("Beginning PS extraction")

powcmbtrue = sum(cmbmap*cmbmap)
cmbnoisemap = zeros((nodecscans, norablocks))
for obs in observations:
    powcmbobs = sum(obs * obs)
    obs = obs * (powcmbtrue/powcmbobs)**0.5
    cmbnoisemap = cmbnoisemap + obs
cmbnoisemap = cmbnoisemap * (1.0/len(observationno))

print(sum(cmbmap*cmbmap))
print(sum(cmbnoisemap*cmbnoisemap))

subplot(1, 2, 1)
imshow(cmbmap)
title("CMB TRUE")
colorbar()

subplot(1, 2, 2)
imshow(cmbnoisemap)
title("CMB 1BIT OBS")
colorbar()

show()

# Plot powerspectrum
mapparams = [512, 512, 2, 2]
k, p, err, h = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap)
k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
k = asarray(k)
p = asarray(p)
err = asarray(err)
err0 = asarray(err0)
k0 = asarray(k0)
p0 = asarray(p0)
p = p*k*(k+1)
p0 = p0*k0*(k0+1)

plot(k0, p0)
xscale("log")
yscale("linear")
xlim((0, 2500))
xlabel(r"Mutlipole moment $l$")
title("Observed Power Spectrum")
ylabel(r"$D_l$")
show()


subplot(1, 2, 1)
plot(k, p, "r", label = "CMB 1Bit")
plot(k0, p0, "b", label = "CMB")
title(r"Powerspectrum")
xlabel(r"Wavevector $k$")
ylabel(r"Power")
xscale("log")
yscale("log")
legend()
subplot(1, 2, 2)
#plot(k0, asarray(p0)-asarray(p), label="Difference")
plot(k0, p/p0, label="Ratio")
title(r"Powerspectrum Ratio")
xlabel(r"Wavevector $k$")
ylabel(r"Ratio")
xscale("linear")
show()