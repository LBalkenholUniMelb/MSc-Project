#--- Make necessary imports

from numpy import *
from pickle import *
from matplotlib.pyplot import *


#--- Read in maps

mapshape = (512, 512)
ranks = 10

cmbnoisemap = zeros(mapshape)
#cmbnoisemap1bit = zeros(mapshape)
#cmbnoisemap2bit = zeros(mapshape)
cmbnoisemap3bit = zeros(mapshape)


for rank in range(ranks):

    cmbnoisemapcurrent = load(open("../../SPT Cloud Code/fixed det noise results 3 bit/8hpp/cmbnoisemap" + str(rank) + ".p", "rb"))
    #cmbnoisemap1bitcurrent = load(open("cmbnoisemap1bit" + str(rank) + ".p", "rb"))
    #cmbnoisemap2bitcurrent = load(open("cmbnoisemap2bitopt" + str(rank) + ".p", "rb"))
    cmbnoisemap3bitcurrent = load(open("../../SPT Cloud Code/fixed det noise results 3 bit/8hpp/cmbnoisemap3bitopt" + str(rank) + ".p", "rb"))

    cmbnoisemap = cmbnoisemap + cmbnoisemapcurrent
    #cmbnoisemap1bit = cmbnoisemap1bit + cmbnoisemap1bitcurrent
    #cmbnoisemap2bit = cmbnoisemap2bit + cmbnoisemap2bitcurrent
    cmbnoisemap3bit = cmbnoisemap3bit + cmbnoisemap3bitcurrent

cmbnoisemap = cmbnoisemap/float(ranks)
#cmbnoisemap1bit = cmbnoisemap1bit/float(ranks)
#cmbnoisemap2bit = cmbnoisemap2bit/float(ranks)
cmbnoisemap3bit = cmbnoisemap3bit/float(ranks)

#--- Save combined maps

dump(cmbnoisemap, open("cmbnoisemap8hppcombined.p", "wb"), -1)
#dump(cmbnoisemap1bit, open("cmbnoisemap1bitcombinedreduced.p", "wb"), -1)
#dump(cmbnoisemap2bit, open("cmbnoisemap2bitcombinedreduced.p", "wb"), -1)
dump(cmbnoisemap3bit, open("cmbnoisemap3bitfixeddetnoise8hppcombined.p", "wb"), -1)
