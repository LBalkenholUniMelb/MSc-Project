#--- IMPORTS

from numpy import *
from matplotlib.pyplot import *
from cosmdigitclasses import *
rc('text', usetex=True)

#--- OBSERVATION PARAMETERS
ranks = range(12)
seeds = range(10)
schemes = ["1bit", "2bithm", "2bitopt"]
pixelnumber = 512
mapparams = [pixelnumber, pixelnumber, 2, 2]
cosm = Cosmologist()


#############################################
#############################################
# CUT SHORT BY READ IN
#############################################
#############################################

# #--- READ IN DATA
# allcmbs = []
# all1bitcmbs = []
# all2bithmcmbs = []
# all2bitoptcmbs = []
#
# for rank in ranks:
#     for seed in seeds:
#
#         cmbfname = "cmbrank" + str(rank) + "seed" + str(seed) + ".txt"
#         cmb1bitfname = "digitisedcmb" + schemes[0] + "nobs20rank" + str(rank) + "seed" + str(seed) + ".txt"
#         cmb2bithmfname = "digitisedcmb" + schemes[1] + "nobs20rank" + str(rank) + "seed" + str(seed) + ".txt"
#         cmb2bitoptfname = "digitisedcmb2bitnobsopt20rank" + str(rank) + "seed" + str(seed) + ".txt"
#
#         cmb = zeros((pixelnumber, pixelnumber))
#         j = 0
#         for row in open(cmbfname):
#             rowvalues = [float(i) for i in row.split()]
#             cmb[j] = rowvalues
#             j += 1
#         allcmbs.append(cmb)
#
#         cmb1bit = zeros((pixelnumber, pixelnumber))
#         j = 0
#         for row in open(cmb1bitfname):
#             rowvalues = [float(i) for i in row.split()]
#             cmb1bit[j] = rowvalues
#             j += 1
#         all1bitcmbs.append(cmb1bit)
#
#         cmb2bithm = zeros((pixelnumber, pixelnumber))
#         j = 0
#         for row in open(cmb2bithmfname):
#             rowvalues = [float(i) for i in row.split()]
#             cmb2bithm[j] = rowvalues
#             j += 1
#         all2bithmcmbs.append(cmb2bithm)
#
#         cmb2bitopt = zeros((pixelnumber, pixelnumber))
#         j = 0
#         for row in open(cmb2bitoptfname):
#             rowvalues = [float(i) for i in row.split()]
#             cmb2bitopt[j] = rowvalues
#             j += 1
#         all2bitoptcmbs.append(cmb2bitopt)
#
#     print("------------------------")
#     print("Completed Rank " + str(rank))
#     print("------------------------")
#
#
# #---- FIND PS RATIO FOR ALL CMBS AND AVERAGE
# all1bitcmbratios = []
# all2bithmcmbratios = []
# all2bitoptcmbratios = []
# k = []
# l = 0
# while l < len(allcmbs):
#
#     cmb = allcmbs[l]
#     cmb1bit = all1bitcmbs[l]
#     cmb2bithm = all2bithmcmbs[l]
#     cmb2bitopt = all2bitoptcmbs[l]
#
#
#     k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmb)
#     k0 = asarray(k0)
#     p0 = asarray(p0)
#     p0 = p0 * k0 * (k0 + 1)
#
#     if l == 0:
#         k = k0
#
#     k1, p1, err1, h1 = cosm.sriniPowerSpectrum(mapparams, cmb1bit)
#     k1 = asarray(k1)
#     p1 = asarray(p1)
#     p1 = p1 * k1 * (k1 + 1)
#     p1 = p1 * ((sum(p0 * p0) / sum(p1 * p1)) ** 0.5)
#     all1bitcmbratios.append(p1/p0)
#
#     k2hm, p2hm, err2hm, h2hm = cosm.sriniPowerSpectrum(mapparams, cmb2bithm)
#     k2hm = asarray(k2hm)
#     p2hm = asarray(p2hm)
#     p2hm = p2hm * k2hm * (k2hm + 1)
#     p2hm = p2hm * ((sum(p0 * p0) / sum(p2hm * p2hm)) ** 0.5)
#     all2bithmcmbratios.append(p2hm/p0)
#
#     k2opt, p2opt, err2opt, h2opt = cosm.sriniPowerSpectrum(mapparams, cmb2bitopt)
#     k2opt = asarray(k2opt)
#     p2opt = asarray(p2opt)
#     p2opt = p2opt * k2opt * (k2opt + 1)
#     p2opt = p2opt * ((sum(p2opt * p2opt) / sum(p2opt * p2opt)) ** 0.5)
#     all2bitoptcmbratios.append(p2opt/p0)
#
#     print("PS DONE FOR: " + str(float(float(l)/float(len(allcmbs)))))
#
#     savetxt("cmbpsind" + str(l) + ".txt", p0)
#     savetxt("cmb1bitpsind" + str(l) + ".txt", p1)
#     savetxt("cmb2bithmpsind" + str(l) + ".txt", p2hm)
#     savetxt("cmb2bitoptpsind" + str(l) + ".txt", p2opt)
#
#
#     l += 1
#
#
# savetxt("allcmbpsk.txt", k)


#############################################
#############################################
# READ IN PS
#############################################
#############################################

k = []
for row in open("allcmbpsk.txt"):
    rowvalues = [float(i) for i in row.split()]
    k.append(rowvalues[0])
k = asarray(k)

allcmbps = []
all1bitcmbratios = []
all2bithmcmbratios = []
all2bitoptcmbratios = []
l = 0
while l < 119:

    cmbps = zeros((1, len(k)))
    j = 0
    for row in open("cmbpsind" + str(l) + ".txt"):
        rowvalues = [float(i) for i in row.split()][0]
        cmbps[0][j] = rowvalues
        j += 1

    cmb1ps = zeros((1, len(k)))
    j = 0
    for row in open("cmb1bitpsind" + str(l) + ".txt"):
        rowvalues = [float(i) for i in row.split()][0]
        cmb1ps[0][j] = rowvalues
        j += 1

    cmb2hmps = zeros((1, len(k)))
    j = 0
    for row in open("cmb2bithmpsind" + str(l) + ".txt"):
        rowvalues = [float(i) for i in row.split()][0]
        cmb2hmps[0][j] = rowvalues
        j += 1

    cmb2optps = zeros((1, len(k)))
    j = 0
    for row in open("cmb2bitoptpsind" + str(l) + ".txt"):
        rowvalues = [float(i) for i in row.split()][0]
        cmb2optps[0][j] = rowvalues
        j += 1
    cmb2optps = cmb2optps * ((sum(cmbps * cmbps) / sum(cmb2optps * cmb2optps)) ** 0.5)


    allcmbps.append(cmbps)
    all1bitcmbratios.append((cmb1ps-cmbps)/cmbps)
    all2bithmcmbratios.append((cmb2hmps-cmbps)/cmbps)
    all2bitoptcmbratios.append((cmb2optps-cmbps)/cmbps)

    l += 1

psavg = zeros(shape(k))
cmb1bitavgratio = zeros(shape(k))
cmb2bithmavgratio = zeros(shape(k))
cmb2bitoptavgratio = zeros(shape(k))

l = 0
while l < 119:

    psavg = psavg + allcmbps[l]
    cmb1bitavgratio = cmb1bitavgratio + all1bitcmbratios[l]
    cmb2bithmavgratio = cmb2bithmavgratio + all2bithmcmbratios[l]
    cmb2bitoptavgratio = cmb2bitoptavgratio + all2bitoptcmbratios[l]

    l += 1

psavg = psavg.flatten() * float(1.0/120.0)
cmb1bitavgratio = cmb1bitavgratio.flatten() * float(1.0/120.0)
cmb2bithmavgratio = cmb2bithmavgratio.flatten() * float(1.0/120.0)
cmb2bitoptavgratio = cmb2bitoptavgratio.flatten() * float(1.0/120.0)

#--- PLOT RESULTS
subplot(1, 2, 1)
plot(k, psavg)
xscale("log")
yscale("log")

subplot(1, 2, 2)
plot(k, cmb1bitavgratio, label="1 Bit")
plot(k, cmb2bithmavgratio, label="2 Bit hm")
plot(k, cmb2bitoptavgratio, label="2 Bit opt")
xscale("log")
ylim([-0.01, 0.01])
legend(loc = 2)


show()