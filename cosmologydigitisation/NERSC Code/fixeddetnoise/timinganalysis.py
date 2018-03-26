from numpy import *

setup = []
obs = []
save = []

fname = "timingresults.txt"
for row in open(fname):
    try:
        tval = float(row)
        if tval < 1:
            save.append(tval)
        elif tval > 1 and tval < 10:
            setup.append(tval)
        elif tval > 10:
            obs.append(tval)
        else:
            print("COULD NOT SORT")
    except:
        pass

setup = max(setup)
obs = max(obs)
save = max(save)

timeest = setup+512.0*obs+save

print("MAX TIME EST (m):")
print(1.2*timeest/60.0/60.0)
print(0.73147271534*60.0)
