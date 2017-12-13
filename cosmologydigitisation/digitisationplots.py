from cosmdigitclasses import *
from digitisationschemes import *
from numpy import *
from matplotlib.pyplot import *
rc('text', usetex=True)
rc("xtick", labelsize = 15)
rc("ytick", labelsize = 15)

signal = normal(0, 1, 11)
t = linspace(0, 10, 11)
pow = signal*signal

digit1 = deepcopy(signal)
digit2hm = deepcopy(signal)
digit2opt = deepcopy(signal)

digitise1bit(digit1)
hmx4 = digitise2bithalfmax(digit2hm)
digitise2bitoptimal(digit2opt, 1)

##########################################
##########################################

plot(t, signal, "k")

title(r"Time-ordered Data", fontsize = 20)
xlabel(r"Time", fontsize = 20)
ylabel(r"Signal", fontsize = 20)

ylim((-3.5, 3.5))
xlim((0, 10))

legend(fontsize = 15)

show()

##########################################
##########################################

plot(t, signal, "k")
plot(t, [0 for i in t], color = "k", linestyle = "--")

title(r"1 Bit Level", fontsize = 20)
xlabel(r"Time", fontsize = 20)
ylabel(r"Signal", fontsize = 20)

ylim((-3.5, 3.5))
xlim((0, 10))


show()

##########################################

plot(t, signal, "k", label = "Signal")
plot(t, digit1, "r", label = "1 Bit")

title(r"1 Bit Digitisation", fontsize = 20)
xlabel(r"Time", fontsize = 20)
ylabel(r"Signal", fontsize = 20)

ylim((-3.5, 3.5))
xlim((0, 10))

legend(fontsize = 15)

show()

##########################################
##########################################

plot(t, signal, "k")
plot(t, [0 for i in t], color = "k", linestyle = "--")
plot(t, [hmx4 for i in t], color = "k", linestyle = "--")
plot(t, [-hmx4 for i in t], color = "k", linestyle = "--")

title(r"2 Bit Half Maximum Levels", fontsize = 20)
xlabel(r"Time", fontsize = 20)
ylabel(r"Signal", fontsize = 20)

ylim((-3.5, 3.5))
xlim((0, 10))


show()

##########################################


plot(t, signal, "k", label = "Signal")
plot(t, digit2hm, "r", label = "2 Bit Half Max")

title(r"2 Bit Half Maximum Digitisation", fontsize = 20)
xlabel(r"Time", fontsize = 20)
ylabel(r"Signal", fontsize = 20)

ylim((-3.5, 3.5))
xlim((0, 10))

legend(fontsize = 15)

show()

##########################################
##########################################

plot(t, signal, "k")
plot(t, [0 for i in t], color = "k", linestyle = "--")
plot(t, [0.9816 for i in t], color = "k", linestyle = "--")
plot(t, [-0.9816 for i in t], color = "k", linestyle = "--")

title(r"2 Bit Optimal Levels", fontsize = 20)
xlabel(r"Time", fontsize = 20)
ylabel(r"Signal", fontsize = 20)

ylim((-3.5, 3.5))
xlim((0, 10))


show()

##########################################

plot(t, signal, "k", label = "Signal")
plot(t, digit2opt, "r", label = "2 Bit Opt")

title(r"2 Bit Optimal Digitisation", fontsize = 20)
xlabel(r"Time", fontsize = 20)
ylabel(r"Signal", fontsize = 20)

ylim((-3.5, 3.5))
xlim((0, 10))

legend(fontsize = 15)

show()