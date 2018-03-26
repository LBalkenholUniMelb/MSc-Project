from numpy import *
from pickle import *
from matplotlib.pyplot import *
from numpy.random import *

testmap = normal(0, 1, (512, 512))
testmap[120:512-120, 120:512-120] = 10


imshow(testmap)
colorbar()
show()

dump(testmap, open("ptestmap.p", "wb"))
dump(testmap, open("ptestmapeff.p", "wb"), -1)

testmap = zeros((512, 512))

imshow(testmap)
colorbar()
show()

testmap = load(open("ptestmapeff.p", "rb"))