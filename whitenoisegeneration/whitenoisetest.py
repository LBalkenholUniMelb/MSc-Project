import whiteNoiseGenerator
from matplotlib.pyplot import *
from numpy import arange, mean

# Observation paramters
obsSeconds = 60
obsFreq = 20

# Create time and noise arrays
t = arange(0, obsSeconds, float(1/obsFreq))
whiteNoise = whiteNoiseGenerator.generatewhitenoise(time=obsSeconds, freq=obsFreq)

# Calculate statistics
av = mean(whiteNoise)

# Plot
plot(t, whiteNoise, "bx", label="Signal")
plot([t[0], t[len(t)-1]], [av, av], "r", label="Mean")
xlabel("Time")
ylabel("Signal")
title("White Noise")
legend()
show()