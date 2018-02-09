# Make necessary imports
from numpy import *
from numpy.random import *
from numpy.fft import fftshift, fftfreq, ifft2
#from numpy import *
from matplotlib.pyplot import *

# Define constants
arcmins2radians = (1.0/60.0)*(pi/180.0)

# Define map parameters
fieldsizearcmins = 4096
pixelsizearcmin = 1
pixelnumber = 1024
df = 1.0/(float(fieldsizearcmins)*arcmins2radians)

# Read in data
filelocation = "plik_plus_r0p01_highell_lensedtotCls.dat"
file = open(filelocation)
l = []
dl = []
for row in file:
    rowvalues = [float(i) for i in row.split()]
    l.append(rowvalues[0])
    dl.append(rowvalues[1])

l = asarray(l)
print(len(l))
dl = asarray(dl)

# Calculate cl
cl = (dl*2.0*pi)/(l*(l+1.0))


# Create 2dCl
cl2d = zeros((pixelnumber, pixelnumber))
pixelspacingrad = float(pixelsizearcmin)*arcmins2radians
kx = 2.0*pi*fftfreq(pixelnumber, pixelspacingrad)
kx2d = tile(kx, (pixelnumber, 1))
ky = 2.0*pi*fftfreq(pixelnumber, pixelspacingrad)
ky2d = transpose(tile(ky, (pixelnumber, 1)))
cl2d = zeros((pixelnumber, pixelnumber))
k2d = sqrt(kx2d*kx2d + ky2d*ky2d)


for y in range(pixelnumber):
    for x in range(pixelnumber):
        k = k2d[y, x]
        ind = argmin(abs(l-k))
        cl2d[y, x] = cl[ind]
    print("PROGRESS:")
    print(float(y)/float(pixelnumber)*100.0)

savetxt("cl2df4096r4.txt", cl2d)


print("COMPLETED CALCULATION")



#------

mapf = sqrt(cl2d)*( normal(0, 1, (pixelnumber, pixelnumber)) + 1.0j*normal(0, 1, (pixelnumber, pixelnumber) ) )

mapr = real( ifft2(mapf)[0:int(pixelnumber/2), 0:int(pixelnumber/2)] )

imshow(mapr, interpolation="none")
colorbar()
show()