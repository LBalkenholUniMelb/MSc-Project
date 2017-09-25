#--- IMPORT

from numpy import *
from matplotlib.pyplot import *
from cosmdigitclasses import *
rc('text', usetex=True)

#--- DEFINE OBSERVATION
pixelnumber = 512

#--- READ IN DATA
# Digitised CMB
schemes = ["1bit", "2bithm", "2bitopt"]
cols = ["k", "y", "r", "b"]
cmbdigit = [zeros((512, 512)) for k in schemes]
l = 0
while l < len(cmbdigit):
    filelocation = "SPT Cloud Code/Schemes Results/digitisedcmb" + schemes[l] + "nobs20.txt"
    file = open(filelocation)
    j = 0
    for row in file:
        rowvalues = [float(i) for i in row.split()]
        cmbdigit[l][j] = rowvalues
        j += 1
    l += 1

# CMB
file = open("CMBMap.txt")
cmbmap = zeros((pixelnumber, pixelnumber))
rowindex = 0
for row in file:
    rowvalues = [float(i) for i in row.split()]
    cmbmap[rowindex] = rowvalues[:pixelnumber]
    rowindex += 1

#--- CALUCLATE STATISTICS
mapparams = [512, 512, 2, 2]
cosm = Cosmologist()
k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
k0 = asarray(k0)
p0 = asarray(p0)
p0 = p0 * k0 * (k0 + 1)
cmbdigitpowerspectra = []
j = 0
while j < len(schemes):
    cmbdigitmap = cmbdigit[j]
    k, p, err, h = cosm.sriniPowerSpectrum(mapparams, cmbdigitmap)
    k = asarray(k)
    p = asarray(p)
    p = p*k*(k+1)
    pownorm = (sum(p0*p0)/sum(p*p))**0.5
    cmbdigitpowerspectra.append(p*pownorm)
    j += 1

show()
close()

#--- PLOT RELEVANT QUANTITIES

#--- IMAGES
schemetitles = [r"1 Bit", r"2 Bit Hm", r"2 Bit Opt"]
fig, axes = subplots(nrows=2, ncols=2)
firstplot = True
p = 0
for ax in axes.flat:
    if firstplot:
        im = ax.imshow(cmbmap)
        ax.set_title("CMB")
        firstplot = False
    else:
        im = ax.imshow(cmbdigit[p-1]-cmbmap, vmin=-1.5*10**(-6), vmax=1.5*10**(-6))
        ax.set_title(schemetitles[p-1] + r" - CMB")
        print(amax(cmbdigit[p-1]-cmbmap))
        print(amin(cmbdigit[p-1]-cmbmap))
    p += 1

# Make an axis for the colorbar on the right side
cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
fig.colorbar(im, cax=cax)
show()

#--- RATIO COMPARISON GENERAL

subplot2grid((1, 2), (0, 0))
plot(k0, p0, cols[0], label = r"CMB")

t = 0
while t < len(schemes):
    plot(k0, cmbdigitpowerspectra[t], cols[t+1], label = schemes[t])
    t += 1

xscale("log")
yscale("log")
xlabel(r"Wavevector $k$")
ylabel(r"$p\times k(k+1)$")
title(r"Powerspectrum")
legend(loc = 3)

subplot2grid((1, 2), (0, 1))
t = 0
while t < len(schemes):
    plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
    t += 1

xscale("log")
yscale("linear")
xlabel(r"Wavevector $k$")
ylabel(r"$P/P_0$")
title(r"Powerspectrum Ratio")
ylim([0.9, 1.1])

show()


#--- RATIO COMPARISON UP CLOSE


ax1 = subplot2grid((3, 3), (0, 0), colspan = 3)
t = 0
while t < len(schemes):
    plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
    t += 1
fill_between([2670, 10000], [0, 10], color = "0.7")
#plot([3*10**3, 3*10**3], [0, 2], lw = 2.0, color = "k")
xscale("log")
yscale("linear")
xlabel(r"Wavevector $k$")
ylabel(r"$P/P_0$")
title(r"Powerspectrum Ratio")
legend(loc = 3)
ylim([0.99, 1.01])
xlim([0, max(k0)])

ax2 = subplot2grid((3, 3), (1, 0), colspan = 3)
t = 0
while t < len(schemes):
    plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
    t += 1
fill_between([5695, 10000], [0, 10], color = "0.7")
xscale("log")
yscale("linear")
xlabel(r"Wavevector $k$")
ylabel(r"$P/P_0$")
title(r"Powerspectrum Ratio starting from $k = 3\times 10^3$")
ylim([0.95, 1.05])
xlim([3*10**3, max(k0)])


ax3 = subplot2grid((3, 3), (2, 0), colspan = 3)
t = 0
while t < len(schemes):
    plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
    t += 1
xscale("log")
yscale("linear")
xlabel(r"Wavevector $k$")
ylabel(r"$P/P_0$")
title(r"Powerspectrum Ratio starting from $k = 6 \times 10^3$")
ylim([0.8, 1.2])
xlim([6*10**3, max(k0)])


show()