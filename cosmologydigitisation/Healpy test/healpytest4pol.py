from numpy import *
from matplotlib.pyplot import *
from healpy import *
from numpy.random import *

def dectotheta(dec):
    return pi*(90.0-dec)/180.0

def ratophi(ra):
    return ra*pi/180.0

# SCAN STRATEGY

NSIDE = 128
scanmap = zeros(nside2npix(NSIDE))

resfield = 10
rastart = 0.0
decstart = 0.0
fieldsize = 600.0 #sq deg

noobservations = 100
readoutfreq = 200.0 #Hz
raspeed = 0.2*float(noobservations)#0.00021/60.0 #deg/s
squarehpp = readoutfreq * sqrt(nside2pixarea(NSIDE, True)) / raspeed
cespoints = int((sqrt(fieldsize)/raspeed)*readoutfreq)
nocesscans = int(sqrt(fieldsize/nside2pixarea(NSIDE, True)))
pixlengthscale = sqrt(nside2pixarea(NSIDE))

phi = ratophi(linspace(rastart, rastart + sqrt(fieldsize), cespoints))
theta = dectotheta(linspace(decstart, decstart + sqrt(fieldsize), nocesscans))

obs = 0
while obs < noobservations:

    obsmod = float(((obs+0.5)/(noobservations/2.0)) - 1.0)

    # recreate scan

    for t in theta:
        for p in phi:
            # add pixel offset
            thit = t + obsmod*(pixlengthscale/2.0)
            phit = p + obsmod*(pixlengthscale/2.0)
            if thit < 0:
                thit += pi
            fieldpixind = ang2pix(NSIDE, thit, phit)
            scanmap[fieldpixind] += 1
    print(obsmod)
    obs += 1

mollview(scanmap)
show()


# CMB TEST

Tcmb = 2.73 #K

Dls_len = np.loadtxt('output_planck_r_0.0_2015_cosmo_lensedCls.dat',usecols=[0,1,2,3,4])
l = np.asarray(Dls_len[:,0])
Dl = np.asarray(Dls_len[:,1])
Cl = ( Tcmb**2. * Dl * 2.0 * np.pi ) / ( l * (l + 1.0) )

l = insert(l, (0, 0), (0, 1))
Dl = insert(Dl, (0, 0), zeros(2))
Cl = insert(Cl, (0, 0), zeros(2))

cmbmap = synfast(Cl, NSIDE)
mollview(cmbmap)
show()

resfield = 1000
rastart = 0.0
decstart = 0.0
fieldsize = 600.0 #sq deg
epsilon = 0.1

scanno = 100
scanstarts = linspace(-1.0*epsilon, epsilon, scanno)
s = 0
while s < scanno:

    rastartadj = rastart * (1.0 + scanstarts[s])
    # another layer for dec...
    phi = ratophi(tile(linspace(rastart, rastart + sqrt(fieldsize), resfield), (resfield, 1)))
    theta = dectotheta(transpose(tile(linspace(decstart, decstart + sqrt(fieldsize), resfield), (resfield, 1))))
    fieldpixind = ang2pix(NSIDE, theta, phi)
    cmbmap[fieldpixind] = UNSEEN

    s += 1


phi = ratophi(tile(linspace(rastart, rastart + sqrt(fieldsize), resfield), (resfield, 1)))
theta = dectotheta(transpose(tile(linspace(decstart, decstart + sqrt(fieldsize), resfield), (resfield, 1))))
fieldpixind = ang2pix(NSIDE, theta, phi)
cmbmap[fieldpixind] = UNSEEN

mollview(cmbmap)
show()


clout = anafast(cmbmap)
lout = arange(len(clout))


print(2.0*((180.0/sqrt(nside2pixarea(NSIDE, True)*60.0))**2.0)/pi)


print(amax(lout))


plot(lout, clout * lout * (lout + 1.0))
plot(l, Dl)
#xscale("log")
xlim((0, amax(lout)))
show()