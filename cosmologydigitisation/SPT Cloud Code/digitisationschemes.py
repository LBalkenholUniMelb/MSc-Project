from numpy import inf, digitize, mean, std, median

#--- 1 Bit

def digitise1bit(signal):
    i = 0
    while i < len(signal):
        if signal[i] > 0:
            signal[i] = 1
        else:
            signal[i] = -1
        i += 1

#--- 2 Bit
def digitise2bithalfmax(signal):
    x1 = -inf
    x2 = - 0.5 * max([abs(max(signal)), abs(min(signal))])
    x3 = 0
    x4 = -x2
    x5 = inf
    outlevels = [-3, -1, 1, 3]
    indices = digitize(signal, [x1, x2, x3, x4, x5])
    for i in range(len(signal)):
        signal[i] = outlevels[indices[i] - 1]

def digitise2bitequalnumbers(signal):
    x1 = -inf
    x5 = inf
    x3 = median(signal)
    upperdata = []
    lowerdata = []
    for i in signal:
        if i >= x3:
            upperdata.append(i)
        else:
            lowerdata.append(i)
    x4 = median(upperdata)
    x2 = median(lowerdata)
    outlevels = [-3, -1, 1, 3]
    indices = digitize(signal, [x1, x2, x3, x4, x5])
    for i in range(len(signal)):
        signal[i] = outlevels[indices[i] - 1]

def digitise2bitoptimal(signal, var):
    x1 = -inf
    x4 = 0.9816 * var
    x2 = -x4
    x3 = 0
    x5 = inf
    outlevels = [-1.51 * var, -0.4528 * var, 0.4528 * var, 1.51 * var]
    indices = digitize(signal, [x1, x2, x3, x4, x5])
    for i in range(len(signal)):
        signal[i] = outlevels[indices[i] - 1]