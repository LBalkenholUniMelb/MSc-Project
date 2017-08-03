from numpy.random import normal

def generatewhitenoise(time = 100, freq = 200, mean = 0, sigma = 1):
    return normal(mean, sigma, time*freq)