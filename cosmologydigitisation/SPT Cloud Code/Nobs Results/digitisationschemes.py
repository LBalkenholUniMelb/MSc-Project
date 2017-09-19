#--- 1 Bit

def digitise1bit(signal):
    i = 0
    while i < len(signal):
        if signal[i] > 0:
            signal[i] = 1
        else:
            signal[i] = -1
        i += 1