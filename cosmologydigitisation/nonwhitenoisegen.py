from scipy.fftpack import *
from scipy import *
from numpy import *
from matplotlib.pyplot import *
from cosmdigitclasses import *

cosm = Cosmologist()

def nonwhitenoise(l):
    lknee = 350.0
    return (1.0 + (l/lknee)**(-8.0/3.0))

fieldsizearcmins = 1024.0
pixelsizearcmin = 2.0
pixelnumber = int(fieldsizearcmins/pixelsizearcmin)
df = 1
mapparams = [512, 512, 2, 2]


pixelspacingrad = (pixelsizearcmin/60.0)*(pi/180.0)
lx = fftfreq(pixelnumber, pixelspacingrad)*2.0*pi
lx2d = tile(lx, (pixelnumber, 1))
ly = fftfreq(pixelnumber, pixelspacingrad)*2.0*pi
ly2d = transpose(tile(ly, (pixelnumber, 1)))
l2d = (lx2d*lx2d+ly2d*ly2d)**0.5
l2d = l2d
l2d = l2d + l2d[1, 0]/2.0
l2d = l2d#[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
imshow(l2d)
colorbar()
show()
# imshow(nwnmsk)
# colorbar()
# show()

nwnmsk = nonwhitenoise(l2d)

pall = 0
kall = 0
realno = 100
i = 0
while i < realno:
    noisefr = normal(0, 1, size=(pixelnumber, pixelnumber))
    noisefi = normal(0, 1, size=(pixelnumber, pixelnumber))
    noisef = noisefr + 1j * noisefi
    noisef = noisef*nwnmsk
    noiser = ifft2(noisef)[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
    noiser = real(noiser)
    k, p, err, h = cosm.sriniPowerSpectrum([int(pixelnumber/2), int(pixelnumber/2), 2, 2], noiser)
    p = asarray(p)

    pall = pall + p
    kall = k

    i += 1

pall = pall/realno
plot(kall, pall)
yscale("log")
show()