from numpy import *
from matplotlib.pyplot import *
from scipy.stats import *
from numpy.random import *
from scipy.special import *

rc('text', usetex=True)
rc("xtick", labelsize = 15)
rc("ytick", labelsize = 15)


hpps = linspace(1, 10, 100)
succlvl = -400
sigmanot = 8.52211548825e-07
successrates = asarray([0.0 for i in hpps])

realno = 1

for real in range(realno):

    i = 0
    for hpp in hpps:
        sigma = 7071

        thcontformula = 0.5 * (1.0 + erf(((succlvl) / (sqrt(2.0) * sigma))))


        #dp = normal(loc = 0, scale = sigma, size = hpp)
        #print(len(dp))
        #succno = sum(asarray(dp < succlvl, dtype = int))
        successrates[i] = thcontformula #successrates[i] + float(succno)/float(hpp)
        i += 1

successrates = successrates/float(realno)
print(successrates)

#subplot(2, 1, 1)
#plot(hpps, [succlvl/(sigmanot*sqrt(j)) for j in hpps], "r")
#ylabel("S/N")
#xscale("log")
#show()

#thcont = norm.cdf(succlvl)
#thcontformula = 0.5*(1.0 + erf(((succlvl)/(sqrt(2.0)*sigma))) )

#subplot(2, 1, 2)
plot(hpps, successrates, "b", label = "Drawn from Distribution")
#plot(hpps, [thcontformula for i in hpps], "r", label = r"Theoretical $\pm 10\%$")
#plot(hpps, [thcontformula+0.1*thcontformula for i in hpps], "r")
#plot(hpps, [thcontformula-0.1*thcontformula for i in hpps], "r")
#ylim((0, 0.05))
ylabel("Success rate", fontsize = 20)
xlabel("Hits per pixel", fontsize = 20)
legend(fontsize = 15)
#xscale("log")
show()


#x = linspace(-5, 5, 100)
#plot(x, norm.cdf(x), "b")
#plot(y, y, "rx")
#xlim((min(x), max(x)))
#ylim((min(x), max(x)))
#show()