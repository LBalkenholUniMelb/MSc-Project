from scipy.fftpack import *
from scipy import *
from numpy import *
from matplotlib.pyplot import *
from cosmdigitclasses import *

cosm = Cosmologist()

def nonwhitenoise(l):
    lknee = 350.0
    return (1.0 + (l/lknee)**(-8.0/3.0))#8/3

fieldsizearcmins = 1024.0
pixelsizearcmin = 2.0
pixelnumber = int(fieldsizearcmins/pixelsizearcmin)
df = 1.0
mapparams = [int(pixelnumber/2), int(pixelnumber/2), pixelsizearcmin, pixelsizearcmin]

pixelspacingrad = (pixelsizearcmin/60.0)*(pi/180.0)
lx = fftfreq(pixelnumber, pixelspacingrad)*2.0*pi
lx2d = tile(lx, (pixelnumber, 1))
ly = fftfreq(pixelnumber, pixelspacingrad)*2.0*pi
ly2d = transpose(tile(ly, (pixelnumber, 1)))
l2d = sqrt(lx2d*lx2d+ly2d*ly2d)
#l2d = l2d
l2d = l2d + l2d[1, 0]/2.0 # hit centre of cl bin to avoid first entry being 0
#l2d = l2d#[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
#imshow(l2d)
#colorbar()
#show()


nwnmsk = nonwhitenoise(l2d)

print(nwnmsk)


factor = sqrt(df/2.0)*nwnmsk

pall = 0
kall = 0
realno = 10
i = 0
while i < realno:
    noisefr = normal(0, 1, size=(pixelnumber, pixelnumber))
    noisefi = normal(0, 1, size=(pixelnumber, pixelnumber))
    noisef = noisefr + 1j * noisefi
    noisef = noisef*factor
    noiser = ifft2(noisef)[0:int(pixelnumber/2), 0:int(pixelnumber/2)]
    noiser = real(noiser)

    #imshow(noiser)
    #colorbar()
    #show()

    k, p, err, h = cosm.sriniPowerSpectrum(mapparams, noiser)
    p = asarray(p)
    p = 2.0 * pixelnumber * pixelnumber * p

    pall = pall + p
    kall = k

    print(i)

    i += 1

pall = pall/realno
plot(kall, pall)
ylabel("Cl")
xlabel("Multipole Moment l")
yscale("linear")
xscale("linear")
ylim((1, 10))
show()