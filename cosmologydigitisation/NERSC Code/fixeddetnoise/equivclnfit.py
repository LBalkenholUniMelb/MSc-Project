#--- Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from cosmdigitclasses import *
from numpy import *
from scipy.signal import convolve2d
from scipy.stats import *
from pickle import *
from scipy.optimize import curve_fit


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
raspeed = 0.0005 #arcmin/s
nodecscans = mappixelnumber
norablocks = mappixelnumber
radatapoints = int(((ralims[1]-ralims[0])/raspeed)*readoutfreq)
compression = int(radatapoints/norablocks)
observationlim = 128
compressions = [800000, 800000, 80000, 8000, 800, 80, 8]
raspeeds = [0.0005, 0.0005, 0.005, 0.05, 0.5, 5.0, 50.0]
obslims = [128, 100, 100, 100, 100, 100, 100]
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

k0 = load(open("results/k.p", "rb"))
p0 = load(open("results/p0.p", "rb"))
p1bit = load(open("results/p1bit800000hpp.p", "rb"))
pn1bit = load(open("results/pn1bit800000hpp.p", "rb"))

#--------------------#
#--------------------#
#--------------------#

pi = 0
lpinchoff = 0
while pi < len(k0):
    if p1bit[pi]/pn1bit[pi] >= 1.02:
        lpinchoff = k0[pi]
        break
    pi += 1


def clnfit(l, cln, A, alpha):
    return cln + ( A/(l**alpha) )

print(lpinchoff)
indpinchoff = argmin(abs(k0 - lpinchoff))

k0end = k0[indpinchoff:]
p0end = p0[indpinchoff:]
pn1bitend = pn1bit[indpinchoff:]
p1bitend = p1bit[indpinchoff:]

popt, pcov = curve_fit(clnfit, k0end, pn1bitend, p0 = [2.99453964e-07, 2.76959427e+25, 8.68086299e+00], maxfev = 10000)
print(popt)
print(sqrt(diag(pcov)))
p1bitendth = asarray([clnfit(j, popt[0], popt[1], popt[2]) for j in k0end])
#plot(k0, p1bit/pn1bit, "k")
#plot(k0end, pn1bitend, "b")
plot(k0end, pn1bitend, "y")
plot(k0end, p1bitendth, "r")
#ylim((0.95, 1.05))

xscale("log")
yscale("log")
show()