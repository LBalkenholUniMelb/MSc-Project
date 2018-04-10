#--- Make necessary imports
from matplotlib.pyplot import *
from numpy import zeros, arange, real, shape, round
from numpy.fft import ifft2, fft2, ifft, fftshift, fftfreq
from cosmdigitclasses import *
from numpy import *
from scipy.signal import convolve2d
from scipy.stats import *
from pickle import *
from scipy.optimize import curve_fit, leastsq

rc('text', usetex=True)
rc("xtick", labelsize = 20)
rc("ytick", labelsize = 20)

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

pixrms = noisedet/sqrt(asarray(obslims, dtype=float) * float(pixelsizearcmin)/asarray(raspeeds)) # muK
thnoisemap = pixrms*float(pixelsizearcmin) # muK arcmin

#--- DEFINE FITS

def clnfit(l, cln, A, alpha):
    return cln + (A / (l ** alpha))

def clnleastsqpenalised(params, x, y):
    penalisation = 0
    if params[0] <= pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians * pixrms[i] * pixrms[i]:
        penalisation = (pixelsizearcmin * pixelsizearcmin * arcmins2radians * arcmins2radians * pixrms[i] * pixrms[i] - params[0]) * 100.0
    return y - clnfit(x, params[0], params[1], params[2]) - penalisation

#--- READ IN MAPS

schemes = ["1bit", "2bit", "3bit"]
k = load(open("results/k.p", "rb"))
p0 = load(open("results/p0.p", "rb"))

clnequivs = empty((2, len(schemes), len(hpp)))

k1bitflatstart = [4000, 3500, 3000, 2500]
k1bitflatstop = 6780
k2bitflatstart = [4000, 3500, 3000, 2500]
k2bitflatstop = 6780

flatstart = 4000
flatstop = 6780
flatstartind = argmin(abs(k - flatstart))
flatstopind = argmin(abs(k - flatstop))

i = 0
while i < len(hpp):

    j = 0
    while j < len(schemes):

        scheme = schemes[j]
        pdigit = load(open("results/p" + str(scheme) + str(hpp[i]) + "hpp.p", "rb"))
        pn = load(open("results/pn" + str(scheme) + str(hpp[i]) + "hpp.p", "rb"))

        if hpp[i] > 800000:

            # Determine Cln via fitting

            pi = 0
            lpinchoff = 0
            while pi < len(k):
                if pdigit[pi] / pn[pi] >= 1.02:
                    lpinchoff = k[pi]
                    break
                pi += 1

            indpinchoff = argmin(abs(k - lpinchoff))
            kend = k[indpinchoff:]
            pnend = pn[indpinchoff:]
            pdigitend = pdigit[indpinchoff:]


            popt, pcov, infodict, mesg, ier = leastsq(clnleastsqpenalised, x0 = [[pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixrms[i]*pixrms[i]*1.05, 2.76959427e+25, 8.68086299e+00]], args = (kend, pdigitend), maxfev=10000, full_output=True)
            #print(pcov)
            #print(poptp[0])
            #print(pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixrms[i]*pixrms[i])
            print(popt[0])
            popt[0] = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixrms[i]*pixrms[i] + mean((pdigitend-pnend)[:argmin(abs(kend - 7200))])
            print(popt[0])
            pcov[0][0] = std((pdigitend-pnend)[:argmin(abs(kend - 7200))])


            #plot(kend, pdigitend-pnend, "k")
            #show()

            #popt, pcov = curve_fit(clnfit, kend, pdigitend, p0=[pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixrms[i]*pixrms[i], 2.76959427e+25, 8.68086299e+00], maxfev=10000)
            #plot(kend, pdigitend, "k")
            #plot(kend, pnend, "b")
            #plot(kend, [pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixrms[i]*pixrms[i] for p in kend], "y", marker = "o")
            #plot(kend, [popt[0] for o in kend], "g")
            #plot(kend, [poptp[0] for m in kend], "m")
            #plot(kend, [clnfit(l, poptp[0], poptp[1], poptp[2]) for l in kend], "m")
            #plot(kend, [clnfit(l, popt[0], popt[1], popt[2]) for l in kend], "g")
            #xscale("log")
            #yscale("log")
            #title(str(hpp[i]) + " // " + str(schemes[j]))
            #ylim((0.98, 1.02))
            #show()

            if sqrt(popt[0])/arcmins2radians - thnoisemap[i] < 0:
                print("HIT")
                clnequivs[0][j][i] = clnfit(amax(k), popt[0], popt[1], popt[2]) #clnequivs[0][j-1][i]#-thnoisecl[i]
            else:
                clnequivs[0][j][i] = popt[0]#-thnoisecl[i]
                clnequivs[1][j][i] = sqrt(pcov[0][0])

        else:

            # Determine Cln via averaging flat tail

            #plot(k, pn, "k")
            #plot(k, pdigit, "r")
            #plot([flatstart, flatstart*1.0001], [10e-10, 10e10], "b")
            #plot([flatstop, flatstop*1.0001], [10e-10, 10e10], "b")
            #xscale("log")
            #yscale("log")
            #show()

            clnequivs[0][j][i] = mean(pdigit[flatstartind:flatstopind])#-thnoisecl[i]
            clnequivs[1][j][i] = std(pdigit[flatstartind:flatstopind])


        j += 1


    i += 1

mapequivs = sqrt(clnequivs)/arcmins2radians



#--- PLOT RESULTS


cols = ["r", "b", "y"]
s = 0
while s < len(schemes):
    plot(hpp, (mapequivs[0][s]-thnoisemap)/thnoisemap, marker = "o", color = cols[s] ,label= r"" + schemes[s])
    s += 1

title(r"Induced Map Noise Level through Digitisation", fontsize = 20)
xlabel(r"Hits per pixel", fontsize = 20)
ylabel(r"$\Delta \sigma_{\mathrm{MAP}} / \sigma_{\mathrm{MAP}} $", fontsize = 20)
xscale("log")
yscale("log")
legend(loc = "upper left", fontsize = 15)
show()