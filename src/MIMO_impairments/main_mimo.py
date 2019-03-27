from symbol_generator import SymbolGenerator
from mimo_channel import MIMOChannel
from test_channel import TestChannel
from equalizer import MIMOEqualizer
from slicer import Slicer
import channel_coefficients

import numpy as np
import matplotlib.pyplot as plt
from plotable import fft_plot

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

class GeneratorConfig:
    n_tx_symbs = 100000
    dualpol = True
    complex_type = True
    os_factor = 4

class ChnConfig:
    brate = 32e9
    taps = 230
    os_factor = 4
    tscale, _ = TestChannel.sinc(taps, brate, os_factor*brate)

    h11 = channel_coefficients.h11
    h12 = channel_coefficients.h12
    h21 = channel_coefficients.h21
    h22 = channel_coefficients.h22

    coefficients = np.array([
        [h11, h12],
        [h21, h22]
    ])
    SNR = None

class EqConfig:
    taps = 256
    mu = 2e-5
    dualpol = True
    os_factor = 4
    coefficients = np.zeros((dualpol+1,2, taps), dtype=np.complex_)
    # coefficients[0,10] = 1


#### TIME MEASUREMENT ####
import timeit
starttime = timeit.default_timer()
#### TIME MEASUREMENT ####

######## INSTANCES #############
gen = SymbolGenerator(GeneratorConfig)
chn = MIMOChannel(ChnConfig)
eq = MIMOEqualizer(EqConfig)
# slicer = Slicer(None)

######## DATAFLOW #############
xn = gen.generate_symbols()
chn_out = chn.process(xn)

e = list()
ye = list()
#adapt(salida canal(entradaequ), simbolos, algun otro kwarg)
eq.adapt(chn_out, xn[:,::gen.os_factor], error=e, ye=ye)
#def adapt(self, data, reference, **kwargs):
#### TIME MEASUREMENT ####
elapsed = timeit.default_timer() - starttime
print('Time elapsed: {0:.2f}'.format(elapsed))
print('Samples per second: {0:.2f}'.format(gen.os_factor*gen.n_tx_symbs/elapsed))
#### TIME MEASUREMENT ####

######## PLOT ###########
# Channel
# Time
# fig, axes = plt.subplots(nrows=2, ncols=4, sharey=True, sharex=True)

# plt.suptitle("Channel")
# for x in range(axes.shape[1]):
#     for y in range(axes.shape[0]):
#         p = x//2
#         title = 'H' + str(y) + str(p) + ' '
#         if x%2 == 0:
#             title += 'Real'
#             axes[y,x].stem(chn.tscale, np.real(chn.coefficients[y,p]))
#         else:
#             title += 'Imaginary'
#             axes[y,x].stem(chn.tscale, np.imag(chn.coefficients[y,p]))
#         axes[y,x].grid()
#         axes[y,x].set_title(title)

# Freq
fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True)

plt.suptitle("Channel Freq Response")
for x in range(axes.shape[1]):
    for y in range(axes.shape[0]):
        title = 'H' + str(y) + str(x)
        axes[y,x].grid()
        axes[y,x].set_title(title)
        fft_plot(chn.coefficients[x, y], chn.os_factor*ChnConfig.brate, ax=axes[y,x])

# # Equalizer
# # Time
# fig, axes = plt.subplots(nrows=2, ncols=4, sharey=True, sharex=True)

# plt.suptitle("Equalizer")
# for x in range(axes.shape[1]):
#     for y in range(axes.shape[0]):
#         p = x//2
#         title = 'H' + str(y) + str(p) + ' '
#         if x%2 == 0:
#             title += 'Real'
#             axes[y,x].stem(np.real(eq.coefficients[y,p]))
#         else:
#             title += 'Imaginary'
#             axes[y,x].stem(np.imag(eq.coefficients[y,p]))
#         axes[y,x].grid()
#         axes[y,x].set_title(title)

# Freq
fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True)

plt.suptitle("Equalizer Freq Response")
for x in range(axes.shape[1]):
    for y in range(axes.shape[0]):
        title = 'H' + str(y) + str(x)
        axes[y,x].grid()
        axes[y,x].set_title(title)
        fft_plot(eq.coefficients[x, y], eq.os_factor*ChnConfig.brate, ax=axes[y,x])

# # Error
# fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, sharex=True)

# plt.suptitle("Error")
# for y in range(axes.shape[0]):
#     title = ''
#     if x%2 == 0:
#         title += 'Real'
#         axes[y].plot([abs(err[0,0])**2 for err in e])
#     else:
#         title += 'Imaginary'
#         axes[y].plot([abs(err[1,0])**2 for err in e])
#     axes[y].grid()
#     axes[y].set_title(title)
# # plt.figure()
# # plt.suptitle("Error H")
# # plt.plot([abs(x[0,0])**2 for x in e])
# # plt.grid()
# # plt.figure()
# # plt.suptitle("Error V")
# # plt.plot([abs(x[1,0])**2 for x in e])
# # plt.grid()

# # Coefficients convergence
# fig, axes = plt.subplots(nrows=2, ncols=4, sharey=True, sharex=True)

# plt.suptitle("Coefficients convergence")
# for x in range(axes.shape[1]):
#     for y in range(axes.shape[0]):
#         p = x//2
#         title = 'H' + str(y) + str(p) + ' '
#         if x%2 == 0:
#             title += 'Real'
#             axes[y,x].plot([np.real(k[y,p]) for k in c])
#         else:
#             title += 'Imaginary'
#             axes[y,x].plot([np.imag(k[y,p]) for k in c])
#         axes[y,x].grid()
#         axes[y,x].set_title(title)
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
fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True)

plt.suptitle("Output convergence")
for x in range(axes.shape[1]):
    for y in range(axes.shape[0]):
        if y == 0:
            title = 'H '
        else:
            title = 'V '
        if x%2 == 0:
            title += 'Real'
            axes[y,x].plot([np.real(k[y]) for k in ye], '*')
        else:
            title += 'Imaginary'
            axes[y,x].plot([np.imag(k[y]) for k in ye], '*')
        axes[y,x].grid()
        axes[y,x].set_title(title)
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

plt.show()
# app = QtGui.QApplication([])

# win = pg.GraphicsWindow(title="PyCommSimulator")
# win.resize(1000,600)
# # win.showFullScreen()
# win.setWindowTitle('PyCommSimulator')

# # Enable antialiasing for prettier plots
# pg.setConfigOptions(antialias=True)

# sample_ws = 100
# sample_ss = 50
# p0 = win.addPlot(title="Output convergence", colspan=2)
# p0.setRange(xRange=(-1,1), yRange=(-1,1), padding=0.5)
# p0.showGrid(x=True, y=True,alpha=0.5)
# p0.setLabel('bottom', 'Real')
# p0.setLabel('left', 'Imaginary')
# curve = p0.plot(pen=None, symbol='o', name="Iteration")
# ptr = 0
# def update():
#     global curve, ptr
#     # curve.setData(np.real(y[100*ptr:100*(ptr+1),0]), np.imag(y[100*ptr:100*(ptr+1),0]), name="asd")
#     curve.setData(np.real(y[sample_ss*ptr:sample_ws+sample_ss*ptr,0]), np.imag(y[sample_ss*ptr:sample_ws+sample_ss*ptr,0]))
#     ptr = ptr + 1 if ptr < (y.shape[0]/sample_ss - sample_ws) else 0
# timer = QtCore.QTimer()
# timer.timeout.connect(update)
# timer.start(10)

# win.nextRow()

# p1 = win.addPlot(title="Coefficients convergence")
# p1.showGrid(x=True, y=True, alpha=0.5)
# p1.setRange(xRange=None, yRange=(-0.04,0.6), padding=0.0)
# p1.setLabel('bottom', 'Iteration')
# p1.setLabel('left', 'Amplitude')
# p1curve = list()
# for i in range(c.shape[1]):
#     p1curve.append(p1.plot(pen=(i,c.shape[1])))
# p1ptr = 0
# def updatep1():
#     global p1curve, p1ptr
#     for i in range(c.shape[1]):
#         p1curve[i].setData(np.real(c[0:sample_ws+sample_ss*p1ptr, i]))

#     p1ptr = p1ptr + 1 if p1ptr < (c.shape[0]/sample_ss - sample_ws) else 0
# timer1 = QtCore.QTimer()
# timer1.timeout.connect(updatep1)
# timer1.start(10)

# p2 = win.addPlot(title="Coefficients convergence")
# p2.showGrid(x=True, y=True, alpha=0.5)
# p2.setRange(xRange=None, yRange=(-0.04,0.6), padding=0.0)
# p2.setLabel('bottom', 'Taps')
# p2.setLabel('left', 'Amplitude')
# p2curve = p2.plot(pen='y', symbol='x')
# p2ptr = 0
# def updatep2():
#     global p2curve, p2ptr
#     p2curve.setData(np.real(c[ptr*sample_ss, :]))

#     p2ptr = p2ptr + 1 if p2ptr < (c.shape[0]/sample_ss - sample_ws) else 0
# timer2 = QtCore.QTimer()
# timer2.timeout.connect(updatep2)
# timer2.start(10)

# if __name__ == '__main__':
#     import sys
#     if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
#         QtGui.QApplication.instance().exec_()
