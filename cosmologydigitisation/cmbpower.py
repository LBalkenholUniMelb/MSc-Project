# Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from scipy.stats import binned_statistic
from numpy.fft import ifft2
from cosmdigitclasses import *
from numpy import *

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

# Transform into map
cmbnoisemap = real(fft.ifft2((realpart + 1j*imagpart)))[0:int(pixelnumber/2), 0:int(pixelnumber/2)] + 2.725
#imshow(cmbnoisemap)
#title("CMB Signal and White Noise")
#colorbar()
#show()

# Feed into Srinis code for powerspectrum
# cosm = Cosmologist()
# cosm.sriniPowerSpectrum([512, 512, 2, 2], cmbnoisemap)
# plot(l, dl, "r")
# xlabel("l")
# ylabel("cl")
# yscale("log")
# xscale("log")
# title("Powerspectrum of CMB and Noise")
# show()


# Recreate Scan Strategy
declims = [0, 1024] #arcmins
ralims = [0, 1024] #arcmins
readoutfreq = 500 #Hz
raspeed = 0.1 #arcmin/s
nodecscans = 512
norablocks = 512
#decstep = declims[1]/nodecscans
#rastep = raspeed/readoutfreq
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
cesscans = zeros((nodecscans, radatapoints))

for d in range(nodecscans):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        # Create relevant tod signals here
        sig = normal(cmbnoisemap[d, ri], 0.0006, compression)
        cesscans[d, rstart:rstop] = sig#digitise1bit(sig, meanvalue = cmbnoisemap[d, ri])
    print(str(d))

# Recompress into map
maprecreated = zeros((nodecscans, norablocks))
for d in range(shape(cesscans)[0]):
    for ri in range(norablocks):
        rstart = ri*compression
        rstop = rstart + compression
        m = mean(cesscans[d, rstart:rstop])
        maprecreated[d, ri] = m
    print(str(d))

# Substract maps
imshow(maprecreated)
title("Recreated CMB and Noise Map")
colorbar()
show()

# Compare powerspectra
mapparams = [512, 512, 2, 2]
kfull, clfull = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap, col = "b")
kcomp, clcomp = cosm.sriniPowerSpectrum(mapparams, maprecreated, col = "None", mark="o")
plot(kfull, clfull, "b", label="Orgininal")
plot(kcomp, clcomp, "r", label="Recreated")
xlabel("l")
ylabel("cl")
xscale("log")
yscale("log")
title("Powerspectrum of CMB and Noise Map, original and recreated")
legend()
show()


# nra=1024
# ndec=1024
# pix_arcmin = 2.0
# pix_size = np.pi/180/60*pix_arcmin
# kx = 2*np.pi*np.fft.fftfreq(ndec,d=pix_size)
# kx2d=np.tile(kx,nra)
# ky = 2*np.pi*np.fft.fftfreq(nra,d=pix_size)
# ky2d=(np.tile(ky,ndec)).T
# k = np.round(np.sqrt(kx2d*kx23+ky2d*ky2d))
# cl2d = cl[k]
#
# factor = np.sqrt(blah *cl2d)
#
# realbit=factor*np.random.normal(size=(nra,ndec))
# imagbit=factor*np.random.normal(size=(nra,ndec))
#
# fft_realization=complex(realbit,imagbit)
# out_map = np.float(np.fft.ifft2(fft_realization))
#
# return out_map[0:ndec/2,0:nra/2]
#




#
# # Turn into 2D power spectrum with noise
# # read in noise
# rabins = 100
# decbins = 100
# df = 1
# gwnfftfile = open("gwnfield/gwn2dfft.txt")
# gwnfft = zeros((rabins, decbins))
# gwncmbfft = zeros((rabins, decbins))
# i = 0
# for row in gwnfftfile:
#     gwnfft[i] = row.split()
#     i += 1
#
# # Plot GWN map fft
# imshow(gwnfft)
# title("GWN Map FFT")
# colorbar()
# show()
#
# # adjust cmb to correct length
# cladjstat = binned_statistic(l, cl, "mean", int(((decbins**2)+(rabins**2))**(0.5)))
# cladj = cladjstat[0]
# ladj = cladjstat[1]
# ladjstep = ladj[1]-ladj[0]
# ladj = [ladj[i]+ladjstep/2 for i in range(len(ladj)-1)]
#
# # Plot compressed cmb power spectrum
# plot(ladj, cladj)
# title("Compressed CMB Powerspectrum")
# xlabel("Multipole moment l")
# ylabel("Cl")
# show()
#
# nra=1024
# ndec=1024
# pix_arcmin = 2.0
# pix_size = np.pi/180/60*pix_arcmin
# kx = 2*np.pi*np.fft.fftfreq(ndec,d=pix_size)
# kx2d=np.tile(kx,nra)
# ky = 2*np.pi*np.fft.fftfreq(nra,d=pix_size)
# ky2d=(np.tile(ky,ndec)).T
# k = np.round(np.sqrt(kx2d*kx23+ky2d*ky2d))
# cl2d = cl[k]
#
# factor = np.sqrt(blah *cl2d)
#
# realbit=factor*np.random.normal(size=(nra,ndec))
# imagbit=factor*np.random.normal(size=(nra,ndec))
#
# fft_realization=complex(realbit,imagbit)
# out_map = np.float(np.fft.ifft2(fft_realization))
#
# return out_map[0:ndec/2,0:nra/2]
#
# # create signal+noise 2d power spectrum
# centre = [int(rabins/2), int(decbins/2)]
# for i
# for y in arange(-int(decbins/2), int((decbins/2)), 1):
#     for x in arange(-int(rabins/2), int((rabins/2)), 1):
#         ind = int(((y**2)+(x**2))**0.5)
#         coeff = cladj[ind]
#         gwncmbfft[centre[0]+y, centre[1]+x] = gwnfft[centre[0]+y, centre[1]+x] * ((coeff*(df**2))/2)**0.5
#
# # plot gwn map and power spectrum
# gwnmap = zeros((rabins, decbins))
# i = 0
# for row in open("gwnfield/gwn2drealspace.txt"):
#     gwnmap[i] = row.split()
#     i += 1
#
# imshow(gwnmap)
# title("GWN Map")
# colorbar()
# show()
#
# cosm.sriniPowerSpectrum([100, 100, 2, 2], gwnmap)
# show()
#
# # Plot CMB*Noise FFT
# imshow(gwncmbfft)
# title("CMB and Noise FFT")
# colorbar()
# show()
#
# # Transform back into real space
# gwncmb = real(ifft2(gwncmbfft))
# # cut out corners
# padding = 5
# gwncmb = gwncmb[5:(shape(gwncmb)[0]-padding), 5:(shape(gwncmb)[1]-padding)]
# imshow(gwncmb)
# colorbar()
# show()
#
# # Feed into Srinis code to get 1D powerspectrum
# mapparams = [90, 90, 1.8, 1.8]
# cosm.sriniPowerSpectrum(mapparams, gwncmb)
# show()
#
# # Compare powerspectra
# # plot(l, cl)
# # xscale("log")
# # show()