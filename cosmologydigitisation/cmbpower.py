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
df = 1.0/(fieldsizearcmins*arcmins2radians)
mapfieldsize = int(fieldsizearcmins/2.0)
mappixelnumber = int(pixelnumber/2.0)
declims = [0, mapfieldsize] #arcmins
ralims = [0, mapfieldsize] #arcmins
readoutfreq = 200.0 #Hz
raspeed = 0.5 #arcmin/s
nodecscans = mappixelnumber
norablocks = mappixelnumber
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
observationlim = 1

# noise level
noisemap = 500.0#3.0 # muK arcmin
noisepix = noisemap/float(pixelsizearcmin)
fsky = (mapfieldsize*mapfieldsize)/(4.0*pi*60.0*60.0*(180.0/pi)**2.0)
noisecl = 4.0*pi*fsky*noisepix*noisepix/float(mappixelnumber*mappixelnumber)

print(noisecl)
gsfses

eta = sqrt(float(compression)) * sqrt(float(observationlim)) * sqrt(float(noisecl) / (float(pixelnumber * pixelnumber)))


def nonwhitenoise(l):
    lknee = 350
    if l == 0:
        return 1.0
    else:
        return (1.0 + (abs(l) / lknee) ** (-8.0 / 3.0))

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
fname = "cl2df" + str(fieldsizearcmins) + "r" + str(pixelsizearcmin) + ".txt"
readincld2 = True
if readincld2:
    y = 0
    for row in open(fname):
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
    savetxt("cl2df1024r1.txt", cl2d)

# imshow(cl2d)
# colorbar()
# show()
#
# TEST WHITE NOISE
#cl2d = ones((pixelnumber, pixelnumber), dtype = float)
#
#print("POWER IN DL:")
#print(sum(dl*conjugate(dl)))
#
#print("POWER IN CL2D:")
#print(sum(cl2d*conjugate(cl2d)))


# Combine noise and cmb element-wise

df = 1.0/(float(fieldsizearcmins)*arcmins2radians)**2.0
factor = sqrt((df/2.0)*cl2d)

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


facadj = 2.0#1.192356 * float(fieldsizearcmins)/2048.0

realno = 1
realind = 0
k0 = 0
k = 0
dltot = 0
dlobstot = 0


while realind < realno:

    re = normal(0, 1, (pixelnumber, pixelnumber))
    im = normal(0, 1, (pixelnumber, pixelnumber))

    realpart = factor * re
    imagpart = factor * im
    cmbfreqspace = (realpart + 1.0j*imagpart)

    # Transform into map
    cmbmap = fft.ifft2(cmbfreqspace)[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
    #print("POWER IN CMB REAL SPACE * 4 * N^2:")
    #print(sum(cmbmap*conjugate(cmbmap))*4.0*pixelnumber*pixelnumber)
    #print(sum(cmbmap*conjugate(cmbmap)))

    cmbmap = real(cmbmap)
    cmbmap = cmbmap * pixelnumber * pixelnumber
    imshow(cmbmap)
    colorbar()
    show()

    # print("BECAUSE WE DISCARD IMAGINARY PART:")
    # print(sum(cmbmap*cmbmap*2.0))
    # print(sum(cmbmap*cmbmap))
    #
    # k, p, perr, h = cosm.sriniPowerSpectrum([int(pixelnumber/2), int(pixelnumber/2), pixelsizearcmin, pixelsizearcmin], cmbmap)
    #
    # xscale("linear")
    # yscale("linear")
    # plot(k, pixelnumber*pixelnumber*p*k*(k+1.0)/(pi))
    # xlabel("Multipole Moment", fontsize = 15)
    # ylabel(r"$D_l [\mu K^2]$", fontsize = 15)
    # xscale("linear")
    # yscale("linear")
    # show()
    #
    #print("ULTIMATE:")
    # powps = sum(p*h)*8.0*float(pixelnumber)*float(pixelnumber)
    #
    # powfac = facadj*float(pixelnumber)*float(pixelnumber)
    #
    #
    # powdiffs.append((powfreq-powps)/powfreq)
    #
    # dl = p*k*(k+1.0)/(2.0*pi)
    # dl = dl*powfac
    # print("Adjusting by:")
    # print(powfac)
    #
    # #dl = dl*4.0e6
    #


# plot(k0, pixelnumber*pixelnumber*p0*k0*(k0+1.0)/(pi), "r")
# plot(linput, dlinput, "b")
# xscale("linear")
# yscale("linear")
# xlim((0, 2500))
# xlabel(r"Mutlipole moment $l$")
# title("Observed Power Spectrum")
# ylabel(r"$D_l$")
# show()

# means = []
# j = 0
# while j < i:
#     t = 0
#     m = 0
#     while t <= j:
#         m += powdiffs[t]
#         t += 1
#     means.append(m/(j+1.0))
#     j += 1
#
# dltot = dltot/noreals
#
# plot(k, dltot)
# xscale("linear")
# yscale("linear")
# title("PS")
# show()
#
# xscale("linear")
# yscale("linear")
# cmbmapavg = cmbmapavg/noreals
# imshow(cmbmapavg)
# colorbar()
# show()
#
# xscale("linear")
# yscale("linear")
# powfreqs = asarray(powfreqs)/2.0
# m = mean(powfreqs)
# s = std(powfreqs)
# plot(range(i), powfreqs, "k")
# plot(range(i), [m for q in range(i)], "b")
# plot(range(i), [m+s for q in range(i)], "b")
# plot(range(i), [m-s for q in range(i)], "b")
# plot(range(i), [sum(factor*factor) for p in range(i)], "r")
# xscale("linear")
# yscale("linear")
# show()
#
#
# plot(range(i), means)
# xscale("linear")
# yscale("linear")
# show()
#
# plot(range(i), powdiffs)
# xscale("linear")
# yscale("linear")
# show()

# plot(k, dltot, "r", label = "Reconstructed")
# plot(linput, dlinput, "b", label = "Input")
# title("Powerspectrum", fontsize = 20)
# xscale("linear")
# yscale("linear")
# xlabel(r"Multipole Moment $l$", fontsize = 20)
# ylabel(r"$l(l+1)C_l/2\pi [\mu \mathrm{K}^2]$", fontsize = 20)
# xlim((0, 2000))
# ylim((0, 6000))
# legend(fontsize = 15)
# show()
#
# subplot(2, 1, 1)
# title("Powerspectrum", fontsize = 25)
# plot(k, dltot, "r", label = "Reconstructed")
# plot(linput, dlinput, "b", label= "Input")
# xlabel(r"Multipole Moment $l$", fontsize = 20)
# ylabel(r"$l(l+1)C_l/2\pi [\mu \mathrm{K}^2]$", fontsize = 20)
# legend(fontsize = 20)
# xscale("linear")
# yscale("linear")
# xlim((0, 2500))
#
# # find where linput and k coincide
# keq = []
# psdiff = []
# for count, wavevec in enumerate(k):
#     ind = argmin(abs(linput-wavevec))
#     keq.append(linput[ind])
#     psdiff.append(dlinput[ind]/dltot[count])
#
# subplot(2, 1, 2)
# plot(keq, psdiff)
# title("Powerspectrum Ratio", fontsize = 25)
# xlabel(r"Multipole Moment $l$", fontsize = 20)
# ylabel(r"Input/Reconstructed", fontsize = 20)
# xscale("linear")
# yscale("linear")
# xlim((0, 2500))
# #print(mean(psdiff))
# #print(std(psdiff))
# show()
#


#print("POWER IN CMB REAL SPACE:")
#print(sum(cmbmap*conjugate(cmbmap))*2)
#print(sum(cmbmap*conjugate(cmbmap)))

#print("POWER IN P:")
#print(hits)
#print(sum(p*hitsc))

#plot(l, dl, "b")

# Recreate Scan Strategy


    observationno = range(observationlim)
    observations = [zeros((nodecscans, norablocks)) for i in observationno]
    cesscans = [zeros((nodecscans, radatapoints)) for i in observationno]

    # Create noise for all observations
    noisel = 2.0*pi*fftfreq(norablocks*compression, arcmins2radians*(pixelsizearcmin)/float(compression))
    noisemask = asarray([nonwhitenoise(q) for q in noisel])

    #noisef = normal(0, 1, norablocks*compression)


    #noise = normal(0, eta, size = (len(observationno), nodecscans, norablocks*compression))

    # Decompose into bolometer signals
    for d in range(nodecscans):


        # create noise for this scan
        #noisef = normal(0, eta**0.5, norablocks*compression)#noisemask*normal(0, 1, norablocks * compression)
        #noisef = noisef#*noisemask
        #noiser = ifft(noisef)
        #noiser = real(noiser)

        noiser = normal(0, eta, norablocks*compression)

        for ri in range(norablocks):
            rstart = ri*compression
            rstop = rstart + compression
            # Add noise
            for obs in observationno:
                tod = cmbmap[d, ri] + noiser[rstart:rstop]

                #todpow = sum(tod*tod)
                # digitise here
                #digitise1bit(tod, 1.0)
                #tod = ((todpow / sum(tod * tod)) ** 0.5) * tod
                #tod = digitise2bithalfmax(tod, 0)


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

    #powcmbtrue = sum(cmbmap*cmbmap)
    cmbnoisemap = zeros((nodecscans, norablocks))
    for obs in observations:
    #    powcmbobs = sum(obs * obs)
    #    obs = obs * (powcmbtrue/powcmbobs)**0.5
        cmbnoisemap = cmbnoisemap + obs
    cmbnoisemap = cmbnoisemap * (1.0/len(observationno))

    cmbnoisemap = cmbnoisemap# + normal(0, 1, shape(cmbnoisemap))

    # DIGITISE MAP

    # for dind in range(mappixelnumber):
    #     for rind in range(mappixelnumber):
    #         if cmbnoisemap[dind, rind] >= 0:
    #             cmbnoisemap[dind, rind] = 1.0
    #         else:
    #             cmbnoisemap[dind, rind] = -1.0
    #
    # imshow(cmbnoisemap)
    # colorbar()
    # show()




    k0, p0, err0, h0 = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin], cmbmap)
    k, p, err, h = cosm.sriniPowerSpectrum([mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin],cmbnoisemap)


    # # Normalise Powerspectrum
    # mmi = 0
    # truepow = 0
    # obspow = 0
    # while mmi < len(k0):
    #     if k0[mmi] >= 350 and k0[mmi] <= 1000:
    #         truepow += p0[mmi]**2.0
    #         obspow += p[mmi]**2.0
    #     mmi += 1
    #
    # p = p * (truepow/obspow)**0.5

    dl = float(pixelnumber) * float(pixelnumber) * p0 * k0 * (k0 + 1.0) / pi
    dltot = dltot + dl
    dlobs = pixelnumber * pixelnumber * p * k * (k + 1.0) / pi
    dlobstot = dlobstot + dlobs


    realind += 1




#
# subplot(1, 2, 1)
# imshow(cmbmap)
# title("CMB TRUE")
# colorbar()
#
# subplot(1, 2, 2)
# imshow(cmbnoisemap)
# title("CMB DIGITISED")
# colorbar()
#
# show()

clobstot= (dlobstot*pi)/(k*(k+1.0))
cltot = (dltot*pi)/(k0*(k0+1.0))

#subplot(1, 2, 1)
#plot(k0, cltot, "b")
plot(k, clobstot, "r")
xscale("log")
yscale("log")
xlabel("Multipole Moment $l$")
ylabel("$C_l$")
title("Powerspectrum")
xlim((0, 3500))
show()

subplot(1, 2, 2)
plot(k0, dlobstot/dltot)
xscale("linear")
yscale("linear")
xlim((0, 2500))
ylim((0.9, 1.1))
show()




# Plot powerspectrum
powfix = 2.0*float(pixelnumber)*float(pixelnumber)
mapparams = [mappixelnumber, mappixelnumber, pixelsizearcmin, pixelsizearcmin]
k, p, err, h = cosm.sriniPowerSpectrum(mapparams, cmbnoisemap)
k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
p = powfix*p*k*(k+1.0)/(2.0*pi)
#p0 = powfix*p0*k0*(k0+1.0)/(2.0*pi)

plot(k0, pixelnumber*pixelnumber*p0*k0*(k0+1.0)/(pi))
xscale("linear")
yscale("log")
xlim((0, 2500))
xlabel(r"Mutlipole moment $l$")
title("Observed Power Spectrum")
ylabel(r"$D_l$")
show()


subplot(1, 2, 1)
plot(k, p, "r", label = "CMB + NOISE")
plot(k0, p0, "b", label = "CMB")
title(r"Powerspectrum")
xlabel(r"Wavevector $k$")
xlim((0, 2500))
ylabel(r"Power")
xscale("linear")
yscale("linear")
legend()
subplot(1, 2, 2)
#plot(k0, asarray(p0)-asarray(p), label="Difference")
plot(k0, p/p0, label="Ratio")
title(r"Powerspectrum Ratio")
xlabel(r"Wavevector $k$")
ylabel(r"Ratio")
xscale("linear")
show()