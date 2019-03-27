from symbol_generator import SymbolGenerator
from test_channel import TestChannel
from equalizer import LMSEqualizer
from slicer import Slicer
import numpy as np
import matplotlib.pyplot as plt
from plotable import fft_plot

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

class GeneratorConfig:
    dualpol = False
    n_tx_symbs = 100000
    complex_type = True

class ChannelConfig:
    taps = 160
    brate = 32e9
    Fs = brate*1.2
    tscale, rcoefficients = TestChannel.sinc(taps, brate, brate*1.1)
    _, icoefficients = TestChannel.sinc(taps, brate, brate*1.2)
    coefficients = np.expand_dims(rcoefficients+1j*icoefficients, axis=0)
    SNR = 15
    # SNR = None

class EqConfig:
    taps = 50
    mu = 2e-4
    dualpol = False
    coefficients = np.zeros((dualpol+1, taps), dtype=np.complex_)
    # coefficients[0,10] = 1


######## DECLARACION #############
gen = SymbolGenerator(GeneratorConfig)
chn = TestChannel(ChannelConfig)
eq = LMSEqualizer(EqConfig)
# slicer = Slicer(None)

######## DATAFLOW #############
xn = gen.generate_symbols()

chn_out = chn.process(xn)

# same convolution with even length arrays introduces a 0 in 1st pos
chn_out = chn_out[:,1:]

e = list()
y = list()
c = list()
eq.adapt(chn_out, xn, error=e, ye=y, coefficients=c)


# Por trabajar en una sola polarizacion:
e = np.array([np.absolute(np.squeeze(x, axis=0))**2 for x in e])


y = np.array([np.squeeze(out, axis=0) for out in y])


c = np.array([np.squeeze(coef, axis=0) for coef in c])

######## PLOT ###########
# Channel
# Time
# plt.figure()
# plt.suptitle("Channel")
# ax = plt.subplot(121)
# plt.title("Real")
# plt.stem(chn.tscale, np.real(chn.coefficients[0]))
# plt.grid()
# plt.subplot(122, sharey=ax)
# plt.title("Imaginary")
# plt.stem(chn.tscale, np.imag(chn.coefficients[0]))
# plt.grid()
# # Freq
# plt.figure()
# plt.suptitle("Channel Freq Response")
# fft_plot(chn.coefficients[0], chn.Fs)
# plt.grid()

# # Equalizer
# # Time
# plt.figure()
# plt.suptitle("Time Domain Equalizer")
# ax = plt.subplot(121)
# plt.title("Real")
# plt.stem(np.real(eq.coefficients[0]))
# plt.grid()
# plt.subplot(122, sharey=ax)
# plt.title("Imaginary")
# plt.stem(np.imag(eq.coefficients[0]))
# plt.grid()
# # Freq
# plt.figure()
# plt.suptitle("Equalizer Freq Response")
# fft_plot(eq.coefficients[0], chn.Fs)
# plt.grid()
# # Error
# plt.figure()
# plt.suptitle("Error")
# plt.plot(e)
# plt.grid()
# # Coefficients convergence
# plt.figure()
# plt.suptitle("Coefficients convergence")
# plt.title("Real")
# plt.plot(np.real(c))
# plt.grid()
# plt.figure()
# plt.title("Imaginary")
# plt.plot(np.imag(c))
# plt.grid()
# Output convergence
# plt.figure()
# plt.suptitle("Output convergence")
# ax = plt.subplot(121)
# plt.title("Real")
# plt.plot(np.real(y), '*')
# plt.grid()
# plt.subplot(122, sharey=ax)
# plt.title("Imaginary")
# plt.plot(np.imag(y), '*')
# plt.grid()
# plt.show()

# plt.figure()
# plt.scatter(np.real(y), np.imag(y), np.arange(len(y)))

app = QtGui.QApplication([])

win = pg.GraphicsWindow(title="PyCommSimulator")
win.resize(1000,600)
# win.showFullScreen()
win.setWindowTitle('PyCommSimulator')

# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)

sample_ws = 100
sample_ss = 50
p0 = win.addPlot(title="Output convergence", colspan=2)
p0.setRange(xRange=(-1,1), yRange=(-1,1), padding=0.5)
p0.showGrid(x=True, y=True,alpha=0.5)
p0.setLabel('bottom', 'Real')
p0.setLabel('left', 'Imaginary')
curve = p0.plot(pen=None, symbol='o', name="Iteration")
ptr = 0
def update():
    global curve, ptr
    # curve.setData(np.real(y[100*ptr:100*(ptr+1),0]), np.imag(y[100*ptr:100*(ptr+1),0]), name="asd")
    curve.setData(np.real(y[sample_ss*ptr:sample_ws+sample_ss*ptr,0]), np.imag(y[sample_ss*ptr:sample_ws+sample_ss*ptr,0]))
    ptr = ptr + 1 if ptr < (y.shape[0]/sample_ss - sample_ws) else 0
timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(10)

win.nextRow()

p1 = win.addPlot(title="Coefficients convergence")
p1.showGrid(x=True, y=True, alpha=0.5)
p1.setRange(xRange=None, yRange=(-0.04,0.6), padding=0.0)
p1.setLabel('bottom', 'Iteration')
p1.setLabel('left', 'Amplitude')
p1curve = list()
for i in range(c.shape[1]):
    p1curve.append(p1.plot(pen=(i,c.shape[1])))
p1ptr = 0
def updatep1():
    global p1curve, p1ptr
    for i in range(c.shape[1]):
        p1curve[i].setData(np.real(c[0:sample_ws+sample_ss*p1ptr, i]))

    p1ptr = p1ptr + 1 if p1ptr < (c.shape[0]/sample_ss - sample_ws) else 0
timer1 = QtCore.QTimer()
timer1.timeout.connect(updatep1)
timer1.start(10)

p2 = win.addPlot(title="Coefficients convergence")
p2.showGrid(x=True, y=True, alpha=0.5)
p2.setRange(xRange=None, yRange=(-0.04,0.6), padding=0.0)
p2.setLabel('bottom', 'Taps')
p2.setLabel('left', 'Amplitude')
p2curve = p2.plot(pen='y', symbol='x')
p2ptr = 0
def updatep2():
    global p2curve, p2ptr
    p2curve.setData(np.real(c[ptr*sample_ss, :]))

    p2ptr = p2ptr + 1 if p2ptr < (c.shape[0]/sample_ss - sample_ws) else 0
timer2 = QtCore.QTimer()
timer2.timeout.connect(updatep2)
timer2.start(10)

if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
