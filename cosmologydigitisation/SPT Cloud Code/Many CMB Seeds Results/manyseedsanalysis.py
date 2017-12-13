#--- IMPORTS

from numpy import *
from matplotlib.pyplot import *
from cosmdigitclasses import *
rc('text', usetex=True)
rc("xtick", labelsize = 15)
rc("ytick", labelsize = 15)

#--- OBSERVATION PARAMETERS
ranks = [i for i in range(24)]
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
# #
# #
# #---- FIND PS RATIO FOR ALL CMBS AND AVERAGE
# all1bitcmbratios = []
# all2bithmcmbratios = []
# all2bitoptcmbratios = []
#
# allps = []
# allps1bit = []
# allps2bithm = []
# allps2bitopt = []
#
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
#     allps.append(p0)
#
#     if l == 0:
#         k = k0
#
#     k1, p1, err1, h1 = cosm.sriniPowerSpectrum(mapparams, cmb1bit)
#     k1 = asarray(k1)
#     p1 = asarray(p1)
#     p1 = p1 * k1 * (k1 + 1)
#     #p1 = p1 * ((sum(p0 * p0) / sum(p1 * p1)) ** 0.5)
#     allps1bit.append(p1)
#     #all1bitcmbratios.append(p1/p0)
#     #
#     k2hm, p2hm, err2hm, h2hm = cosm.sriniPowerSpectrum(mapparams, cmb2bithm)
#     k2hm = asarray(k2hm)
#     p2hm = asarray(p2hm)
#     p2hm = p2hm * k2hm * (k2hm + 1)
#     #p2hm = p2hm * ((sum(p0 * p0) / sum(p2hm * p2hm)) ** 0.5)
#     allps2bithm.append(p2hm)
#     #all2bithmcmbratios.append(p2hm/p0)
#
#     k2opt, p2opt, err2opt, h2opt = cosm.sriniPowerSpectrum(mapparams, cmb2bitopt)
#     k2opt = asarray(k2opt)
#     p2opt = asarray(p2opt)
#     p2opt = p2opt * k2opt * (k2opt + 1)
#     #p2opt = p2opt * ((sum(p0 * p0) / sum(p2opt * p2opt)) ** 0.5)
#     allps2bitopt.append(p2opt)
#     #all2bitoptcmbratios.append(p2opt/p0)
#
#     print("PS DONE FOR: " + str(float(float(l)/float(len(allcmbs)))))
# #
#     savetxt("cmbpsind" + str(l) + ".txt", p0)
#     savetxt("cmb1bitpsind" + str(l) + ".txt", p1)
#     savetxt("cmb2bithmpsind" + str(l) + ".txt", p2hm)
#     savetxt("cmb2bitoptpsind" + str(l) + ".txt", p2opt)
# #
# #
#     l += 1
# #
# #
# savetxt("allcmbpsk.txt", k)
#
#
# ps = zeros(shape(allps[0]))
# ps1bit = zeros(shape(allps1bit[0]))
# ps2bithm = zeros(shape(allps1bit[0]))
# ps2bitopt = zeros(shape(allps1bit[0]))
# i = 0
# while i < len(allps):
#     ps = ps + allps[i]
#     ps1bit = ps1bit + allps1bit[i]
#     ps2bithm = ps2bithm + allps2bithm[i]
#     ps2bitopt = ps2bitopt + allps2bitopt[i]
#     i += 1
#
# ps = ps * 1.0/float(len(allps))
# ps1bit = ps1bit * 1.0/float(len(allps1bit))
# ps2bitopt = ps2bitopt * 1.0/float(len(allps2bithm))
# ps2bithm = ps2bithm * 1.0/float(len(allps2bitopt))
#
# plot(k, ps, label="True")
# plot(k, ps1bit, label = "1")
# plot(k, ps2bithm, label="2HM")
# plot(k, ps2bitopt, label="2OPT")
# legend()
# xscale("log")
# yscale("log")
# show()


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
allcmbps1 = []
allcmbps2hm = []
allcmbps2opt = []
all1bitcmbratios = []
all2bithmcmbratios = []
all2bitoptcmbratios = []
all1bitcmbdiffs = []
all2bithmcmbdiffs = []
all2bitoptcmbdiffs = []
l = 0
while l < 240:

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
    #cmb2optps = cmb2optps * ((sum(cmbps * cmbps) / sum(cmb2optps * cmb2optps)) ** 0.5)


    allcmbps.append(cmbps)
    allcmbps1.append(cmb1ps)
    allcmbps2hm.append(cmb2hmps)
    allcmbps2opt.append(cmb2optps)
    #all1bitcmbratios.append(cmb1ps)#/cmbps)
    #all2bithmcmbratios.append(cmb2hmps)#/cmbps)
    #all2bitoptcmbratios.append(cmb2optps)#/cmbps)
    #all1bitcmbdiffs.append((cmb1ps - cmbps))
    #all2bithmcmbdiffs.append((cmb2hmps - cmbps))
    #all2bitoptcmbdiffs.append((cmb2optps - cmbps))

    l += 1

psavg = zeros(shape(k))
ps1avg = zeros(shape(k))
ps2hmavg = zeros(shape(k))
ps2optavg = zeros(shape(k))
# cmb1bitavgratio = zeros(shape(k))
# cmb2bithmavgratio = zeros(shape(k))
# cmb2bitoptavgratio = zeros(shape(k))
# cmb1bitavgdiff = zeros(shape(k))
# cmb2bithmavgdiff = zeros(shape(k))
# cmb2bitoptavgdiff = zeros(shape(k))

l = 0
while l < 240:

    psavg = psavg + allcmbps[l]

    ps1avg = ps1avg + allcmbps1[l]
    ps2hmavg = ps2hmavg + allcmbps2hm[l]
    ps2optavg = ps2optavg + allcmbps2opt[l]

    # cmb1bitavgratio = cmb1bitavgratio + all1bitcmbratios[l]
    # cmb2bithmavgratio = cmb2bithmavgratio + all2bithmcmbratios[l]
    # cmb2bitoptavgratio = cmb2bitoptavgratio + all2bitoptcmbratios[l]
    # cmb1bitavgdiff = cmb1bitavgdiff + all1bitcmbdiffs[l]
    # cmb2bithmavgdiff = cmb2bithmavgdiff + all2bithmcmbdiffs[l]
    # cmb2bitoptavgdiff = cmb2bitoptavgdiff + all2bitoptcmbdiffs[l]

    l += 1

psavg = psavg.flatten() * float(1.0/240.0)
ps1avg = ps1avg.flatten() * float(1.0/240.0)
ps2hmavg = ps2hmavg.flatten() * float(1.0/240.0)
ps2optavg = ps2optavg.flatten() * float(1.0/240.0)

# cmb1bitavgratio = cmb1bitavgratio.flatten() * float(1.0/240.0)
# cmb2bithmavgratio = cmb2bithmavgratio.flatten() * float(1.0/240.0)
# cmb2bitoptavgratio = cmb2bitoptavgratio.flatten() * float(1.0/240.0)
# cmb1bitavgdiff = cmb1bitavgdiff.flatten() * float(1.0/240.0)
# cmb2bithmavgdiff = cmb2bithmavgdiff.flatten() * float(1.0/240.0)
# cmb2bitoptavgdiff = cmb2bitoptavgdiff.flatten() * float(1.0/240.0)



#############################################
#############################################
# STARTING SCHEMES 1 SEED IMPORT
#############################################
#############################################
#
# #--- READ IN DATA
# # Digitised CMB
# schemes = ["1bit", "2bithm", "2bitopt"]
# cols = ["k", "y", "r", "b"]
# cmbdigit = [zeros((512, 512)) for k in schemes]
# l = 0
# while l < len(cmbdigit):
#     filelocation = "../Schemes Results/digitisedcmb" + schemes[l] + "nobs20.txt"
#     file = open(filelocation)
#     j = 0
#     for row in file:
#         rowvalues = [float(i) for i in row.split()]
#         cmbdigit[l][j] = rowvalues
#         j += 1
#     l += 1
#
# # CMB
# file = open("../../CMBMap.txt")
# cmbmap = zeros((pixelnumber, pixelnumber))
# rowindex = 0
# for row in file:
#     rowvalues = [float(i) for i in row.split()]
#     cmbmap[rowindex] = rowvalues[:pixelnumber]
#     rowindex += 1
#
# #--- CALUCLATE STATISTICS
# mapparams = [512, 512, 2, 2]
# cosm = Cosmologist()
# k0, p0, err0, h0 = cosm.sriniPowerSpectrum(mapparams, cmbmap)
# k0 = asarray(k0)
# p0 = asarray(p0)
# p0 = p0 * k0 * (k0 + 1)
# cmbdigitpowerspectra = []
# j = 0
# while j < len(schemes):
#     cmbdigitmap = cmbdigit[j]
#     k, p, err, h = cosm.sriniPowerSpectrum(mapparams, cmbdigitmap)
#     k = asarray(k)
#     p = asarray(p)
#     p = p*k*(k+1)
#     pownorm = (sum(p0*p0)/sum(p*p))**0.5
#     cmbdigitpowerspectra.append(p*pownorm)
#     j += 1
#
# show()
# close()
#
#--- PLOT RELEVANT QUANTITIES

#--- IMAGES
schemetitles = [r"1 Bit", r"2 Bit Hm", r"2 Bit Opt"]
cols = ["k", "y", "r", "b"]
#--- RATIO COMPARISON GENERAL

# subplot2grid((1, 2), (0, 0))
# plot(k0, p0, cols[0], label = r"CMB")
#
# t = 0
# while t < len(schemes):
#     plot(k0, cmbdigitpowerspectra[t], cols[t+1], label = schemes[t])
#     t += 1
#
# xscale("log")
# yscale("log")
# xlabel(r"Wavevector $k$")
# ylabel(r"$p\times k(k+1)$")
# title(r"Powerspectrum")
# legend(loc = 3)
#
# subplot2grid((1, 2), (0, 1))
# t = 0
# while t < len(schemes):
#     plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
#     t += 1
#
# xscale("log")
# yscale("linear")
# xlabel(r"Wavevector $k$")
# ylabel(r"$P/P_0$")
# title(r"Powerspectrum Ratio")
# ylim([0.9, 1.1])
#
# show()


#--- RATIO COMPARISON UP CLOSE


# ax1 = subplot2grid((3, 3), (0, 0), colspan = 3)
# t = 0
# while t < len(schemes):
#     plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
#     t += 1
# fill_between([2670, 10000], [0, 10], color = "0.7")
# #plot([3*10**3, 3*10**3], [0, 2], lw = 2.0, color = "k")
# xscale("log")
# yscale("linear")
# xlabel(r"Wavevector $k$")
# ylabel(r"$P/P_0$")
# title(r"Powerspectrum Ratio")
# legend(loc = 3)
# ylim([0.99, 1.01])
# xlim([0, max(k0)])
#
# ax2 = subplot2grid((3, 3), (1, 0), colspan = 3)
# t = 0
# while t < len(schemes):
#     plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
#     t += 1
# fill_between([5695, 10000], [0, 10], color = "0.7")
# xscale("log")
# yscale("linear")
# xlabel(r"Wavevector $k$")
# ylabel(r"$P/P_0$")
# title(r"Powerspectrum Ratio starting from $k = 3\times 10^3$")
# ylim([0.95, 1.05])
# xlim([3*10**3, max(k0)])
#
#
# ax3 = subplot2grid((3, 3), (2, 0), colspan = 3)
# t = 0
# while t < len(schemes):
#     plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
#     t += 1
# xscale("log")
# yscale("linear")
# xlabel(r"Wavevector $k$")
# ylabel(r"$P/P_0$")
# title(r"Powerspectrum Ratio starting from $k = 6 \times 10^3$")
# ylim([0.8, 1.2])
# xlim([6*10**3, max(k0)])
#
#
# show()


#############################################
#############################################
# FINISHING SCHEMES ONE SEED IMPORT
#############################################
#############################################
#
# #--- PLOT 1 SEED VS MANY SEEDS RATIO
#
# subplot2grid((1, 2), (0, 0))
# t = 0
# while t < len(schemes):
#     plot(k0, cmbdigitpowerspectra[t]/p0, cols[t+1], label = schemes[t])
#     t += 1
#
# xscale("log")
# yscale("log")
# xlabel(r"Multipole Moment $l$", fontsize = 20)
# ylabel(r"Observed Power / True Power", fontsize = 20)
# title(r"Powerspectrum Ratio of True and Observed CMB", fontsize = 20)
# #ylim([0.9, 1.1])
# legend(loc=2)
#
# subplot2grid((1, 2), (0, 1))
# plot(k, cmb1bitavgratio, "y", label="1 Bit")
# plot(k, cmb2bithmavgratio, "r", label="2 Bit hm")
# plot(k, cmb2bitoptavgratio, "b", label="2 Bit opt")
#
# xscale("log")
# yscale("log")
# xlabel(r"Mutlipole Moment $l$", fontsize = 20)
# ylabel(r"Observed Power / True Power", fontsize = 20)
# title(r"Powerspectrum Ratio for 240 CMB Seeds", fontsize = 20)
# #ylim([0.9, 1.1])
# legend(loc=2)
# show()
#
#
# #--- PLOT RESULTS
# plot(k, psavg)
# xscale("log")
# yscale("log")
# title(r"Average Power Spectrum")
# xlabel(r"Wavevector $k$")
# ylabel(r"$p \times k (k + 1)$")
# show()
#
# # subplot2grid((2, 2), (0, 0), colspan = 2)
# # plot(k, cmb1bitavgdiff, label="1 Bit")
# # plot(k, cmb2bithmavgdiff, label="2 Bit hm")
# # plot(k, cmb2bitoptavgdiff, label="2 Bit opt")
# # title(r"Avergae Difference")
# # xlabel(r"Wavevector $k$")
# # xscale("log")
# # #ylim([-0.01, 0.01])
# # legend(loc = 2)
#
# subplot2grid((1, 2), (0, 0), colspan = 2)
# fill_between([10, 10**5], [0.999, 0.999], [1.001, 1.001], color="0.7", label=r"$\pm 0.1\%$")
#
# plot(k, cmb1bitavgratio, "y", label="1 Bit")
# plot(k, cmb2bithmavgratio, "r", label="2 Bit hm")
# plot(k, cmb2bitoptavgratio, "b", label="2 Bit opt")
# legend(loc=2)
# title(r"Average Ratio")
# xlabel(r"Wavevector $k$")
# ylim([0.99, 1.01])
# xlim([10, 10**4])
# xscale("log")
#
# show()

#############################################
#############################################
# ANALYSIS
#############################################
#############################################

subplot(1, 2, 1)

plot(k, psavg, "k", label = "CMB")
plot(k, ps1avg, "y", label = "1 Bit")
plot(k, ps2hmavg, "r", label = "2 Bit HM")
plot(k, ps2optavg, "b", label = "2 Bit Opt")
legend(loc = 3)
xscale("log")
yscale("log")


subplot(1, 2, 2)

plot(k, ps1avg/psavg, "y", label = "1 Bit")
plot(k, ps2hmavg/psavg, "r", label = "2 Bit HM")
plot(k, ps2optavg/psavg, "b", label = "2 Bit Opt")


show()