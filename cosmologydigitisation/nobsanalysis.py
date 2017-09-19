#--- IMPORT

from numpy import *
from matplotlib.pyplot import *
from cosmdigitclasses import *
rc('text', usetex=True)

#--- DEFINE OBSERVATION
pixelnumber = 512

#--- READ IN DATA
# Digitised CMB
allnobs = [1] + list(arange(10, 210, 10))
cmbdigit = [zeros((512, 512)) for k in allnobs]
l = 0
for nobs in allnobs:
    filelocation = "SPT Cloud Code/Nobs Results/digitisedcmbnobs" + str(nobs) + ".txt"
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
while j < len(allnobs):
    cmbdigitmap = cmbdigit[j]
    print(j)
    k, p, err, h = cosm.sriniPowerSpectrum(mapparams, cmbdigitmap)
    k = asarray(k)
    p = asarray(p)
    p = p*k*(k+1)
    pownorm = (sum(p0*p0)/sum(p*p))**0.5
    cmbdigitpowerspectra.append(p*pownorm)
    j += 1

show()
close()

t = 0
meandiffs = []
firstnobgr1 = True
while t < len(allnobs):
    meandiffs.append(mean(cmbdigitpowerspectra[t]-p0))
    lbl = r"$N_{obs} > 1$"
    col = "r"
    if t == 0:
        lbl = "$N_{obs} = 1$"
        col = "b"
        plot(k0, cmbdigitpowerspectra[t] - p0, col, label=lbl)
    else:
        if firstnobgr1:
            plot(k0, cmbdigitpowerspectra[t]-p0, col, label=lbl)
            firstnobgr1 = False
        else:
            plot(k0, cmbdigitpowerspectra[t]-p0, col)

    t += 1

#plot(allnobs, meandiffs)
title(r"PS Difference btwn CMB and Digitised Observations", fontsize = 20)
xlabel(r"Multipole Moment $k$", fontsize = 15)
ylabel(r"$\Delta P \times k(k+1)$", fontsize = 15)
xscale("log")
#yscale("log")
legend()
show()