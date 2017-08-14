########################################
########################################
# NECESSARY IMPORTS
########################################
########################################

from numpy.random import normal
from matplotlib.pyplot import *
from numpy import arange, mean, any, zeros, add, nan_to_num, digitize, argwhere, isnan, ones, array_equal, shape, inf, median, std, conj, real, linspace, amax, amin
from scipy.stats import binned_statistic_2d
from copy import deepcopy
from numpy.linalg import norm
from numpy.fft import fft, fftfreq, fftshift, fft2
from pow_spec import *

########################################
########################################
# SIGNAL
# Carries a given signal for a time, readout frequency and a short description
########################################
########################################

class Signal:

    # Initialiser
    ########################################
    def __init__(self, time, freq, signal, description = "Signal", meanvalue = 0, var = 0):
        self.time = time
        self.freq = freq
        self.signal = signal
        self.description = description

        # If mean and variance are npt supplied calculate them here
        if meanvalue == 0 and var == 0:
            self.meanvalue = mean(self.signal)
            self.var = std(self.signal)
        else:
            self.meanvalue = meanvalue
            self.var = var

    # Plot the signal
    ########################################
    def plot(self, showmean = False, showvar = 0, colours = [], showlegend = False, showplot = True):
        # Need to be precise about colours supplied
        # Check what colour arguments are supplied and fill if necessary
        if len(colours) == 0:
            colours = ["k"]
            if showmean:
                colours.append("r")
            for i in range(showvar):
                colours.append("m")
        elif (len(colours) == 1 and (showmean == True or showvar != 0)):
            colours.append("r")
            for i in range(showvar):
                colours.append("m")

        # Plot according to arguments supplied
        t = arange(0, self.time, float(1 / self.freq))
        plot(t, self.signal, colours[0], label=self.description)
        av = self.meanvalue
        if showmean:
            plot([t[0], t[len(t) - 1]], [av, av], colours[1], label="Mean")
        if showvar != 0:
            var = self.var
            for k in range(showvar):
                plot([t[0], t[len(t) - 1]], [av+(k+1)*var, av+(k+1)*var], colours[2 + k], label=str(k + 1) + " STD")
                plot([t[0], t[len(t) - 1]], [av-(k+1)*var, av-(k+1)*var], colours[2 + k])

        xlabel("Time")
        ylabel("Signal")
        title(self.description)
        if showlegend:
            legend()
        if showplot:
            show()

    # Old normalisation methods
    def normalise(self, tc = 1):
        coeff = self.normalisationCoefficient()
        if coeff != 0:
            self.signal = [i * (tc/coeff) for i in self.signal]

    def normalisationCoefficient(self):
        return norm(self.signal)


########################################
# GWN Subclass
########################################

class GaussianWhiteNoise(Signal):

    # Initialiser
    ########################################
    def __init__(self, time, freq, meanvalue, var):
        gwnsignal = normal(meanvalue, var, time * freq)
        Signal.__init__(self, time, freq, gwnsignal, "GWN", meanvalue=meanvalue, var=var)

########################################
########################################
# OBSERVATION
# Carries CES data and bolometer signal
########################################
########################################

class Observation:

    # Calculate RA
    def calculatera(self):
        # Do a proper calculation based on telescope postition, observation time here
        return self.phi

    # Calculate DEC
    def calculatedec(self):
        # Do a proper calculation based on telescope postition, observation time here
        return self.theta

    # Initialiser
    ########################################
    def __init__(self, starttime, endtime, readoutfrequency, phi, theta, bolometersignal):
        self.starttime = starttime
        self.endtime = endtime
        self.readoutfrequency = readoutfrequency
        self.phi = phi
        self.theta = theta
        self.bolometersignal = bolometersignal
        self.ra = self.calculatera()
        self.dec = self.calculatedec()

    # Plot Observation data
    ########################################
    def plotObservation(self):
        # Plot key observation parameters
        time = arange(self.starttime, self.endtime, float(1/self.readoutfrequency))
        plot(time, self.phi, label="phi")
        plot(time, self.theta, label="theta")
        xlabel("Time")
        ylabel("Phi/Theta")
        title("Observation Data")
        legend()
        show()

########################################
########################################
# MAP
# A CMB map, initialised via an array of observations spanning a path of the sky
########################################
########################################
class Map:

    # Check dimension of the area covered by observations
    ########################################
    def checkDimensions(self, observations):
        ramin = observations[0].ra[0]
        ramax = ramin
        decmin = observations[0].dec[0]
        decmax = decmin

        for obs in observations:
            if min(obs.ra) < ramin:
                ramin = min(obs.ra)
            elif max(obs.ra) > ramax:
                ramax = max(obs.ra)

            if min(obs.dec) < decmin:
                decmin = min(obs.dec)
            elif max(obs.dec) > decmax:
                decmax = max(obs.dec)

        return [[ramin, ramax], [decmin, decmax]]

    # Calculates bin edges
    ########################################
    def calculateBinEdges(self, rares, decres):
        # Including the right edge of the last bin
        rastep = float((self.dimensions[0][1] - self.dimensions[0][0]) / rares)
        decstep = float((self.dimensions[1][1] - self.dimensions[1][0]) / decres)
        rabins = [self.dimensions[0][0] + i * rastep for i in range(rares+1)]
        decbins = [self.dimensions[1][0] + i * decstep for i in range(decres+1)]
        return [rabins, decbins]

    # Initialiser
    ########################################
    def __init__(self, observations, rares, decres):
        # Assign instance variables
        self.observations = observations
        self.rares = rares
        self.decres = decres
        self.dimensions = self.checkDimensions(self.observations)
        self.map = zeros((self.rares, self.decres))
        self.bins = self.calculateBinEdges(self.rares, self.decres)

        obsbinnednotnorm = []

        # Add time averaged bolometer signals into the correct bin
        for obs in self.observations:
            obsbinned, xedges, yedges, binno = binned_statistic_2d(obs.ra, obs.dec, obs.bolometersignal.signal, "mean", [self.bins[0], self.bins[1]])
            obsbinnednotnorm.append(obsbinned)
            self.map = self.map + nan_to_num(obsbinned)

        # Normalise to account for different observations falling into partially overlapping bins

        # Check if arrays have the same nan position; NEED TO OPTIMISE THIS
        obsinbin = ones((self.rares, self.decres), dtype=int)
        for i in range(len(obsbinnednotnorm)):
            for j in range(len(obsbinnednotnorm)):
                if i != j:
                    if sameNan(obsbinnednotnorm[i], obsbinnednotnorm[j]):
                        nonNanPositions = findNonNanPositions(obsbinnednotnorm[i])
                        for nonNanPos in nonNanPositions:
                            obsinbin[nonNanPos[0], nonNanPos[1]] += 1

        # Normalise values; NEED TO OPTIMISE THIS
        for i in range(self.rares):
            for j in range(self.decres):
                n = obsinbin[i, j] # total hits
                if n != 1:
                    n = n - 1 # excess hits
                    t = 0.5+(0.25+n)**(0.5) # adequate multiplication value from n=t(t+1)
                    self.map[i, j] = self.map[i, j] * (1/t)

    # Plot the map
    ########################################
    def plotMap(self, noise = False, intp = None):
        imshow(self.map, interpolation = intp)
        titlestring = "CMB Temperature Map"
        if noise:
            titlestring = "White Noise Map"
        title(titlestring)
        xlabel("RA Pixel Number")
        ylabel("DEC Pixel Number")
        colorbar()
        show()


########################################
########################################
# DIGITISER
# Digitises a given signal according to different schemes, returns new signal instance, leaving the old signal unchanged
########################################
########################################
class Digitiser:

    # Initialiser
    ########################################
    def __init__(self):
        pass

    # 1 Bit Digitisation
    ########################################
    def digitise1bit(self, signalin, meanvalue = 0):
        signalout = deepcopy(signalin)
        for i in range(len(signalout.signal)):
            if signalout.signal[i] >= meanvalue:
                signalout.signal[i] = meanvalue+1
            else:
                signalout.signal[i] = meanvalue-1
        signalout.description = signalout.description + " 1Bit"
        return signalout

    # 2 Bit Digitisation, Different schemes available, can return threshold values
    ########################################
    def digitise2bit(self, signalin, scheme, meanvalue = 0, var = 1, xlevels = False):
        # Digitise Signal into 2 Bit, returns new signal
        # Operate on zeromean signal
        signalout = deepcopy(signalin)
        zeromeansignal = [i - meanvalue for i in signalout.signal]
        des = ""
        x1 = -inf
        x2 = 0
        x3 = 0
        x4 = 0
        x5 = inf
        outlevels = [-3, -1, 1, 3]

        # Half maximum scheme
        if scheme == "halfmax":
            x4 = 0.5*max([abs(max(zeromeansignal)), abs(min(zeromeansignal))])
            x3 = 0
            x2 = -x4
            des = " 2Bit Hm"

        # Roughly equal numbers in each bin
        elif scheme == "equalnumbers":
            x3 = median(zeromeansignal)
            upperdata = []
            lowerdata = []
            for i in zeromeansignal:
                if i >= x3:
                    upperdata.append(i)
                else:
                    lowerdata.append(i)
            x4 = median(upperdata)
            x2 = median(lowerdata)
            des = " 2Bit Eq"

        # Optimal scheme according to Max (1960)
        elif scheme == "optimal":
            x4 = 0.9816*var
            x2 = -x4
            outlevels = [-1.51*var, -0.4528*var, 0.4528*var, 1.51*var]
            des = " 2Bit Opt"

        # Shift according to mean and digitise given the levels supplied by the specific scheme
        bins = [x1, x2, x3, x4, x5]
        indices = digitize(zeromeansignal, bins)
        for i in range(len(zeromeansignal)):
            signalout.signal[i] = outlevels[indices[i] - 1]+ meanvalue
        signalout.description = signalout.description + des
        if xlevels:
            return signalout, [i + meanvalue for i in bins]
        else:
            return signalout


########################################
########################################
# COSMOLOGIST
# Calculates various CMB statistics, compares digitised Signals
########################################
########################################
class Cosmologist:

    # Initialiser
    ########################################
    def __init__(self):
        pass

    # Compare Signal to 1 Bit Digit
    ########################################
    def compareGwn1Bit(self, gwn, showmean = False, showvar = 0, colours = [[], []]):
        # Need to be precise about colour specification
        gwn.plot(showmean = showmean, showvar = showvar, colours = colours[0], showlegend = False, showplot = False)
        signaldigit = Digitiser()
        gwn1bit = signaldigit.digitise1bit(gwn, gwn.meanvalue)
        gwn1bit.plot(showmean = False, showvar = 0, colours = colours[1], showlegend = False, showplot = False)
        title("GWN Signal and 1Bit Digitisation")
        legend()
        show()

    # Compare Signal to 2 Bit Digit
    ########################################
    def compareGwn2BitDigit(self, gwn, scheme = "all", showmean = False, showvar = 0, colours = [[], [], [], []]):
        # Can compare single schemes or all, need to be precise about colour specification
        gwn.plot(showmean=showmean, showvar=showvar, colours=colours[0], showlegend=False, showplot=False)
        signaldigit = Digitiser()
        if scheme == "all":
            gwn2bitopt = signaldigit.digitise2bit(gwn, "optimal", gwn.meanvalue, gwn.var)
            gwn2bithalf = signaldigit.digitise2bit(gwn, "halfmax", gwn.meanvalue, gwn.var)
            gwn2biteq = signaldigit.digitise2bit(gwn, "equalnumbers", gwn.meanvalue, gwn.var)
            gwn2bitopt.plot(showmean=False, showvar=0, colours=colours[1], showlegend=False, showplot=False)
            gwn2bithalf.plot(showmean=False, showvar=0, colours=colours[2], showlegend=False, showplot=False)
            gwn2biteq.plot(showmean=False, showvar=0, colours=colours[3], showlegend=False, showplot=False)
            title("GWN Signal and 2Bit Digitisations")
            legend()
            show()
        else:
            gwn2bitdigit = signaldigit.digitise2bit(gwn, scheme, gwn.meanvalue, gwn.var)
            gwn2bitdigit.plot(showmean=False, showvar=0, colours=colours[1], showlegend=False, showplot=False)
            title("GWN Signal and 2Bit " + scheme + " Digitisation")
            legend()
            show()

    # Old Power Spectrum codes
    ########################################
    # def calculateBolometerSignalPowerSpectrum(self, signal, showmean = False, showvar = False, onlypositivel = True):
    #     # Calculate the Powerspectrum of a Bolometer Signal
    #     sig = signal.signal
    #     time = signal.time
    #     freq = signal.freq
    #     spacing = 1/freq
    #     datapoints = time*freq
    #     powerspec = [real(i*conj(i)) for i in fft(sig)]
    #     l = fftfreq(datapoints, spacing)
    #     lpos = []
    #     powerspecpos = []
    #     if onlypositivel:
    #         for i in range(len(l)):
    #             if l[i] > 0:
    #                 lpos.append(l[i])
    #                 powerspecpos.append(powerspec[i])
    #     else:
    #         lpos = l[1:]
    #         powerspecpos = powerspec[1:]
    #     semilogy(lpos, powerspecpos, "kx", label = "Powerspectrum")
    #     showlegend = False
    #     if showmean or showvar:
    #         meanvalue = mean(powerspecpos)
    #         showlegend = True
    #         if showmean:
    #             plot(lpos, [meanvalue for i in lpos], "r", label = "Mean")
    #         if showvar:
    #             var = std(powerspecpos)
    #             plot(lpos, [meanvalue + var for i in lpos], "m", label = "Mean +- Std")
    #             plot(lpos, [meanvalue - var for i in lpos], "m")
    #     if showlegend:
    #         legend()
    #     title("Powerspectrum")
    #     xlabel("Wavenumber l")
    #     ylabel("mod al sq")
    #     show()
    #
    # def calculateBolometerSignalCrossPowerSpectrum(self, signals, showmean = False, showvar = False, onlypositivel = True):
    #     # Calculate the Cross Power of Two Bolometer Signals
    #     signal1 = signals[0]
    #     signal2 = signals[1]
    #     sig1 = signal1.signal
    #     sig2 = signal2.signal
    #     time = signal1.time
    #     freq = signal1.freq
    #     spacing = 1/freq
    #     datapoints = time*freq
    #     crosspowerspec = [real(fft(sig1)[k]*conj(fft(sig2)[k])) for k in range(datapoints)]
    #     l = fftfreq(datapoints, spacing)
    #     lpos= []
    #     crosspowerspecpos = []
    #     if onlypositivel:
    #         for i in range(len(l)):
    #             if l[i] > 0:
    #                 lpos.append(l[i])
    #                 crosspowerspecpos.append(crosspowerspec[i])
    #     else:
    #         lpos = l[1:]
    #         crosspowerspecpos = crosspowerspec[1:]
    #     plot(lpos, crosspowerspecpos, "kx", "Crosspower")
    #     showlegend = False
    #     if showmean or showvar:
    #         showlegend += True
    #         meanvalue = mean(crosspowerspecpos)
    #         if showmean:
    #             plot(lpos, [meanvalue for i in lpos], "r", label = "Mean")
    #         if showvar:
    #             var = std(crosspowerspecpos)
    #             plot(lpos, [meanvalue + var for i in lpos], "m", label = "Mean +- Std")
    #             plot(lpos, [meanvalue - var for i in lpos], "m")
    #     if showlegend:
    #         legend()
    #     title("Crosspowerspectrum")
    #     xlabel("Wavenumber l")
    #     ylabel("Crosspower")
    #     show()
    #
    # def calculate2DPowerSpectrum(self, map):
    #     # Calculate the 2D Power Spectrum of a Map
    #     powerspec = fft2(map.map)
    #     powerspec = fftshift(powerspec)
    #     powerspec = powerspec * conj(powerspec)
    #     powerspec = powerspec.real
    #     imshow(powerspec, norm=LogNorm())
    #     colorbar()
    #     show()
    #
    # def calculate1DPowerSpectrum(self, map1, map2 = None, maskingfunction = None, maskingfunctionparameters = None, showmean = False, showvar = False):
    #     # Calculates the 1D Powerspectrum of an Image
    #     # 2D
    #     if maskingfunction == None:
    #         def maskingfunction(r, p):
    #             return 1
    #     crosspower = True
    #     if map2 == None:
    #         crosspower = False
    #         map2 = map1
    #     powerspec1 = fft2(map1.map)
    #     powerspec1 = fftshift(powerspec1)
    #     powerspec2 = fft2(map2.map)
    #     powerspec2 = fftshift(powerspec2)
    #     powerspec = powerspec1 * conj(powerspec2)
    #     powerspec = powerspec.real
    #     # Azimuthal Average
    #     ymax, xmax = shape(powerspec)
    #     center = [xmax/2.0, ymax/2.0]
    #     r = zeros(shape(powerspec))
    #     # NEED TO OPTIMISE THIS
    #     for x in range(xmax):
    #         for y in range(ymax):
    #             r[y, x] = ((x-center[0])**2 + (y-center[1])**2)**(0.5)
    #             powerspec[y, x] = powerspec[y, x]*maskingfunction(r[y, x], maskingfunctionparameters)
    #     bins = linspace(amin(r), amax(r), int(xmax/2))
    #     indices = digitize(r, bins)
    #     powerspec1D = [0 for i in range(amax(indices))]
    #     powerspec1Derr = [0 for i in range(amax(indices))]
    #     powerspec1Dcount = [0 for i in range(amax(indices))]
    #     for x in range(xmax):
    #         for y in range(ymax):
    #             powerspec1D[indices[y, x]-1] += powerspec[y, x]
    #             powerspec1Dcount[indices[y, x]-1] += 1
    #     for i in range(len(powerspec1D)):
    #         if powerspec1Dcount[i] != 0:
    #             powerspec1D[i] = powerspec1D[i]/powerspec1Dcount[i]
    #             powerspec1Derr[i] = powerspec1D[i]/((powerspec1Dcount[i])**(0.5))
    #     # Cut zero freq
    #     powerspec1D = powerspec1D[1:]
    #     powerspec1Derr = powerspec1Derr[1:]
    #     lbins = range(len(powerspec1D))
    #     # Plot
    #     errorbar(lbins, powerspec1D, yerr = powerspec1Derr)
    #     showlegend = False
    #     if showmean or showvar:
    #         showlegend = True
    #         meanvalue = mean(powerspec1D)
    #         if showmean:
    #             plot(lbins, [meanvalue for i in lbins], "r", label="Mean")
    #         if showvar:
    #             var = std(powerspec1D)
    #             plot(lbins, [meanvalue + var for i in lbins], "m", label="Mean +- Std")
    #             plot(lbins, [meanvalue - var for i in lbins], "m")
    #     if showlegend:
    #         legend()
    #     title("Crosspowerspectrum")
    #     if not crosspower:
    #         title("Powerspectrum")
    #     xlabel("Wavenumberbin")
    #     ylabel("Power")
    #     show()

    # Srinis code to get (Cross) Powerspectrum
    ########################################
    def sriniPowerSpectrum(self, mapparams, map1, map2 = None, col = "b", mark = "None"):
        use_white_noise = 0  # 0 : CMB; 1 - white noise
        clf()
        if use_white_noise:
            if map2 == None:
                ax = subplot(111, yscale='log')  # do not use logx here
            else:
                ax = subplot(111)
        else:
            ax = subplot(111, xscale='log', yscale='log')
        if not use_white_noise:
            #CMB_UNLEN = pickle.load(gzip.open('CMB_lensed.pkl.gz', 'rb'))  # unlensed CMB
            #tSZ = pickle.load(gzip.open('tSZ.pkl.gz', 'rb'))  # tSZ of a 2e15 cluster with beta model
            #CMB_UNLEN_with_tSZ = pickle.load(gzip.open('CMB_lensed.pkl.gz', 'rb'))  # sum of the unlensed CMB and beta model
            #CMB_LEN = pickle.load(gzip.open('CMB_unlensed.pkl.gz','rb'))  # lensed CMB by a 2e15 cluster with NFW profile (slightly smoothed acoustic peaks but difficult to see for 1 cluster)
            #CMB_DIFF = CMB_LEN - CMB_UNLEN  # difference between unlensed and lensed CMB
            #noofsims = 10  # len(CMB_LEN)
            noofsims = 1
        else:
            noofsims = 1

        for n in range(noofsims):
            #print('Sim = %s' % (n + 1))
            if not use_white_noise:
                # read the CMB map
                # IMAGE = np.loadtxt('CMB_sim_0.txt') #CMB

                # use any of the below CMB and compute the power
                # IMAGE = CMB_UNLEN[n]
                # IMAGE = CMB_LEN[n]
                # IMAGE = CMB_DIFF[n]
                # IMAGE = CMB_UNLEN_with_tSZ[n]
                # clf();imshow(CMB_DIFF[n]);colorbar();show();sys.exit()

                IMAGE = map1

                # these are the dimension of all CMB sims
                nx, ny = IMAGE.shape
                dx = dy = 2.  # arcmins
            else:
                # similar to what prof. Lopez asked. 1000 elements; 6 x 6 arcsec^2 box
                # nx = ny = 1000;dx = dy = 0.0001 #arcminutes

                nx, ny, dx, dy = mapparams
                IMAGE = np.random.standard_normal((nx, ny))
                IMAGE_CROSS = None
                if map2 != None:
                    IMAGE_CROSS_L = map2.map
                    IMAGE_CROSS_L -= np.mean(IMAGE_CROSS_L)
                    IMAGE_CROSS = np.copy(IMAGE_CROSS_L)
                if n == noofsims - 1:
                    IMAGE_L = map1.map
                    IMAGE_L -= np.mean(IMAGE_L)
                    IMAGE = np.copy(IMAGE_L)
                    # nx, ny = IMAGE.shape
                    # dx, dy = 2.,2.

                    # clf();imshow(IMAGE);colorbar();show();sys.exit()

            #IMAGE -= np.mean(IMAGE)  # data has zero mean like CMB
            #clf();imshow(IMAGE, cmap = cm.jet,interpolation='None');colorbar();show();sys.exit();#savefig('CMB.png', bbox_inches='tight', pad_inches=0.);sys.exit()

            mapparams = [nx, ny, dx, dy]
            powspec = 0
            if map2 == None:
                powspec = fn_plot_pow_spec(mapparams, IMAGE)  # auto power spectrum of the image
            else:
                powspec = fn_plot_pow_spec(mapparams, IMAGE, IMAGE_CROSS)  # cross power spectrum of the images
            k, p_k, p_k_err, hitcount = zip(*powspec)

            # powspec = fn_plot_pow_spec(mapparams, MAP1 = IMAGE, MAP2 = tSZ) #crosspower between CMB_UNLEN_with_tSZ and tSZ
            # k, p_k, p_k_err = zip(*powspec); p_k = np.abs(p_k)

            # p_k/=max(p_k) #normalising to unity
            if n == 0:
                avg_p_k = np.copy(p_k)
            elif n < noofsims - 1:
                avg_p_k += np.copy(p_k)

            # errorbar(k, p_k, yerr = [p_k_err, p_k_err], marker = '.', color = 'gray', alpha = 0.1)#, label = 'Sim = %s' %(n+1))
            if n == noofsims - 1:
                lb = "Powerspectrum"
                if map2 != None:
                    lb = "Crosspowerspectrum"
                cl = [p_k[i]*k[i]*(k[i]+1)/(2*np.pi) for i in range(len(k))]
                #plot(k[1:], cl[1:], marker=mark, color=col, label=lb)  # , label = 'Sim = %s' %(n+1))
                #av = mean(p_k)
                #var = std(p_k)
                #plot(k, [av for i in k], marker='None', color='k', label='Mean')
                #plot(k, [av+var for i in k], marker='None', color='gray', label='Mean +- Std')
                #plot(k, [av-var for i in k], marker='None', color='gray')

            else:
                plot(k, p_k, marker='None', color='gray', alpha=0.5)  # , label = 'Sim = %s' %(n+1))

        #avg_p_k /= (noofsims - 1)
        #plot(k, avg_p_k, 'r.', label='Average')
        #legend(loc=1, fancybox=1)  # ;ylim(1e-5,1.)
        if map2 == None:
            title("Powerspectrum")
        else:
            title("Crosspowerspectrum")
        #xlabel('Multipole k');
        #ylabel('Power')
        return k, p_k, p_k_err, hitcount
        #legend(fancybox = 1)
        #show()


########################################
########################################
# HELPER FUNCTIONS
########################################
########################################
# Determines coinciding Nan positions
def sameNan(A, B):
    # Check if arrays have the same NaN positions
    return array_equal(argwhere(isnan(A)), argwhere(isnan(B)))

# Returns indicies of non Nan positions
def findNonNanPositions(A):
    # return an array indicating the non-Nan positions
    nanPos = argwhere(isnan(A))
    nonNanArray = ones(shape(A), dtype = bool)
    for pos in nanPos:
        nonNanArray[pos[0], pos[1]] = False
    return argwhere(nonNanArray)

# Cutoff function, first parameter is cutoff value
def cutOff(r, p):
    if r < p[0]:
        return 1
    else:
        return 0

# Digitissation schemes, change the input signal
def digitise1bit(signal, meanvalue = None):
    if meanvalue == None:
        meanvalue = mean(signal)
    for i in signal:
        if i >= meanvalue:
            i = meanvalue+1
        else:
            i = meanvalue-1
    return signal

def digitise2bithalfmax(signal, meanvalue = None):
    if meanvalue == None:
        meanvalue = mean(signal)
    signal = signal - meanvalue
    x1 = -inf
    x2 = - 0.5 * max([abs(max(signal)), abs(min(signal))])
    x3 = 0
    x4 = -x2
    x5 = inf
    outlevels = [-3, -1, 1, 3]
    indices = digitize(signal, [x1, x2, x3, x4, x5])
    for i in range(len(signal)):
        signal[i] = outlevels[indices[i] - 1] + meanvalue

def digitise2bitequalnumbers(signal, meanvalue = None):
    if meanvalue == None:
        meanvalue = mean(signal)
    signal = signal - meanvalue
    x1 = -inf
    x5 = inf
    x3 = median(signal)
    upperdata = []
    lowerdata = []
    for i in signal:
        if i >= x3:
            upperdata.append(i)
        else:
            lowerdata.append(i)
    x4 = median(upperdata)
    x2 = median(lowerdata)
    outlevels = [-3, -1, 1, 3]
    indices = digitize(signal, [x1, x2, x3, x4, x5])
    for i in range(len(signal)):
        signal[i] = outlevels[indices[i] - 1] + meanvalue

def digitise2bitoptimal(signal, meanvalue = None, var = None):
    if meanvalue == None:
        meanvalue = mean(signal)
    if var == None:
        var = std(signal)
    x1 = -inf
    x4 = 0.9816 * var
    x2 = -x4
    x3 = 0
    x5 = inf
    outlevels = [-1.51 * var, -0.4528 * var, 0.4528 * var, 1.51 * var]
    indices = digitize(signal, [x1, x2, x3, x4, x5])
    for i in range(len(signal)):
        signal[i] = outlevels[indices[i] - 1] + meanvalue
