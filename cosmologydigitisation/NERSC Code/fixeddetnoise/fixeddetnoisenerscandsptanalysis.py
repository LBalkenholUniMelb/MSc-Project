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

sizeofskysqarcmin = 4.0*(pi*(180.0/pi)**2.0)*60.0*60.0
fsky = float(mapfieldsize*mapfieldsize)/sizeofskysqarcmin

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

#p0 = ones(shape(k))
#p0 = asarray([i for i in k])


clnequivs = empty((3, len(schemes), len(hpp)))

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
        k = load(open("results/k.p", "rb"))

        # --- REBIN K

        deltak = k[1] - k[0]
        newk = linspace(amin(k) - deltak / 2.0, amax(k) + deltak / 2.0, (len(k) / 10) + 1)

        newkinds = digitize(k, newk[:len(newk) - 1])

        newpn = zeros(shape(newk[:len(newk) - 1]))
        newpdigit = zeros(shape(newk[:len(newk) - 1]))
        hits = histogram(newkinds, len(newk[:len(newk) - 1]))[0]

        print(schemes[j])
        print(hpp[i])
        if raw_input("Display PS ratio? (y/)") == "y":
            plot(k, pdigit/pn, "k")
            #plot(k, pdigit, "r")
            #yscale("log")
            #xscale("log")
            show()
            close()

        for counter, newind in enumerate(newkinds):
            newpdigit[newind - 1] += pdigit[counter]
            newpn[newind - 1] += pn[counter]


        for counter, h in enumerate(hits):
            newpdigit[counter] = newpdigit[counter] / float(h)
            newpn[counter] = newpn[counter] / float(h)

        pn = newpn
        pdigit = newpdigit

        k = newk[:len(newk) - 1] + ((newk[1] - newk[0]) / 2.0)
        flatstart = 4000
        flatstop = 6780
        flatstartind = argmin(abs(k - flatstart))
        flatstopind = argmin(abs(k - flatstop))


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

            clne = pixelsizearcmin*pixelsizearcmin*arcmins2radians*arcmins2radians*pixrms[i]*pixrms[i] + mean((pdigitend-pnend))
            deltaclne = 2.0*std((pdigitend-pnend)[:argmin(abs(kend - 7200))])/sqrt(float(len((pdigitend-pnend))))

            clnequivs[0][j][i] = clne
            clnequivs[1][j][i] = deltaclne
            #clnequivs[2][j][i] = sqrt(2.0/((2.0*flatstart+1)*fsky*float(len((pdigitend-pnend)[:argmin(abs(kend - 7200))]))))*popt[0]


        else:

            # Determine Cln via averaging flat tail

            clnequivs[0][j][i] = mean(pdigit[flatstartind:])
            clnequivs[1][j][i] = 2.0*std(pdigit[flatstartind:])/sqrt(float(len(pdigit[flatstartind:])))
            #clnequivs[2][j][i] = sqrt(2.0 / ((2.0 * flatstart + 1) * fsky * float(len((pdigitend - pnend)[:argmin(abs(kend - 7200))])))) * clnequivs[0][j][i]

        j += 1

    i += 1

mapequivs = sqrt(clnequivs)/arcmins2radians
mapequivs[1] = mapequivs[0]*0.5*(clnequivs[1]/clnequivs[0])
mapequivs[2] = mapequivs[0]*0.5*(clnequivs[2]/clnequivs[0])

#--- PLOT RESULTS


cols = ["r", "b", "y"]
s = 0
while s < len(schemes):
    errorbar(hpp, (mapequivs[0][s]-thnoisemap)/thnoisemap, yerr = ((1.0/thnoisemap)-(1.0/mapequivs[0][s]))*mapequivs[1][s], marker = "o", color = cols[s] ,label= r"" + schemes[s])
    s += 1

title(r"Induced Map Noise Level through Digitisation", fontsize = 20)
xlabel(r"Hits per pixel", fontsize = 20)
ylabel(r"$\Delta \sigma_{\mathrm{MAP}} / \sigma_{\mathrm{MAP}} $", fontsize = 20)
xscale("log")
yscale("log")
legend(loc = "upper left", fontsize = 15)
ylim((0.01, 1))
show()