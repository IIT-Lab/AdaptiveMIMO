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
from adaptive_LMS_equalizer import LMS_equalizer as LMSequ
from polarization_generator import Polarization
from DSP_utils import *

def plot_constellations(noise_on, pol, r1_quad, r2_conj_quad, r1_skew, r2_conj_skew, r1_impairments, r2_conj_impairments):
    if noise_on:
        str_title = '+ noise'
    else:
        str_title = 'without noise'

    plt.figure()
    plt.title('Rx constellation' + pol + ' : IQ imbalance ' + str_title)
    plt.scatter(np.real(r1_quad), np.imag(r2_conj_quad), marker='.', c='lime')
    plt.grid(linestyle='--')
    plt.ylim(-2, 2)
    plt.xlim(-2, 2)

    plt.figure()
    plt.title('Rx constellation'  + pol +' : Skew ' + str_title)
    plt.scatter(np.real(r1_skew), np.imag(r2_conj_skew), marker='.', c='fuchsia')
    plt.grid(linestyle='--')
    plt.ylim(-2, 2)
    plt.xlim(-2, 2)

    plt.figure()
    plt.title('Rx constellation'  + pol + 'skew + IQ imbalance ' + str_title)
    plt.scatter(np.real(r1_impairments), np.imag(r2_conj_impairments), marker='.', c='r')
    plt.grid(linestyle='--')
    plt.ylim(-2, 2)
    plt.xlim(-2, 2)

    plt.show()






SNR_to_simulate = 9
noise = False
filter_taps = 15






'''Genero primer polarizacion, i.e, dos canales complejos correspondientes a polarizacion X'''
X_polarization = Polarization()
X_polarization.generate_lanes(noise,SNR_to_simulate)
'''Genero segunda polarizacion, i.e, dos canales complejos correspondientes a polarizacion Y'''
Y_polarization = Polarization()
Y_polarization.generate_lanes(noise,SNR_to_simulate)

'''Notacion: 
    - dos vectores de simbolos ak,i, con i = {1,2}, para cada portadora, en polarizacion X
    - dos vectores de simbolos bk,i, con i = {1,2}, para cada portadora, en polarizacion Y 

'''

#Notacion:     senal_subcarrier_polarizacion_ADICIONAL...
# Ej  r2_y_conj_skew:    senal compleja recibida correspondiente al segundo canal (o segunda subcarrier)
#                        en polarizacion y, y conjugada, con el efecto de skew incluido
#Vectores de simbolos complejos generados
ak_x_ch1 = X_polarization.symb_ch1_complex
bk_y_ch1 = Y_polarization.symb_ch1_complex

ak_x_ch2 = X_polarization.symb_ch2_complex
bk_y_ch2 = Y_polarization.symb_ch2_complex

#Vectores con unicamente skew incluido
r1_x_skew = X_polarization.r1_skew
r1_y_skew = Y_polarization.r1_skew

r2_x_conj_skew = X_polarization.r2_conj_skew
r2_y_conj_skew = Y_polarization.r2_conj_skew

#Vectores con IQ imbalance incluido
r1_x_quad = X_polarization.r1_quad
r1_y_quad = Y_polarization.r1_quad

r2_x_conj_quad = X_polarization.r2_conj_quad
r2_y_conj_quad= Y_polarization.r2_conj_quad

#Vectores con todos los impairments incluidos
r1_x_impairments = X_polarization.r1_impairments
r1_y_impairments = Y_polarization.r1_impairments

r2_x_conj_impairments = X_polarization.r2_conj_impairments
r2_y_conj_impairments = Y_polarization.r2_conj_impairments


equalizer_taps = 15
len_filter = np.copy(equalizer_taps)
#step = 0.0001
step = 0.0045
delay_equalizer = int((len_filter-1)/2)
# Creamos instancia del ecualizador
LMS_equalizer = LMSequ(equalizer_taps, step)


'''Pruebo con primer canal, polarizacion x'''
#Obtengo simbolos enviados
real_ak = np.real(ak_x_ch1)
imag_ak = np.imag(ak_x_ch1)

real_signal_imp = np.real(r1_x_impairments)
imag_signal_imp = np.imag(r1_x_impairments)


out_equ_x, err_x = LMS_equalizer.adapt(real_signal_imp)

real_resh = real_ak[delay_equalizer:(len(real_ak)-(delay_equalizer+1))]
#real_resh = np.sign(real_resh)


detected = [0] * len(real_ak)
#symbols_ref = np.array([[-1. + 1j * 1.], [-1. - 1j * 1.], [1. + 1j * 1.0], [1.0 - 1j * 1.0]])
symbols_ref = np.array([[-1], [1]])

for muestra in range(len(out_equ_x)):
    argmin = np.argmin((np.abs(out_equ_x[muestra] - symbols_ref)) ** 2)
    detected[muestra] = symbols_ref[argmin, 0]

count = 0
for i in range (len(real_resh)):
    if (detected[i]!= real_ak[i]):
        count +=1

print('Para ' +str(X_polarization.n_data) + 'datos')
print('Hubo ' +str(count) + ' errores')


# plt.figure()
# plt.scatter(real_ak, np.zeros(len(real_ak)), marker='.', c='lime')
# plt.grid(linestyle='--')
# plt.ylim(-2, 2)
# plt.xlim(-2, 2)
# plt.show()
#
# plt.figure()
# plt.scatter(real_signal_imp, np.zeros(len(real_ak)), marker='.')
# plt.grid(linestyle='--')
# plt.ylim(-2, 2)
# plt.xlim(-2, 2)
# plt.show()
#
# plt.figure()
# plt.scatter(out_equ_x, np.zeros(len(out_equ_x)), marker='.')
# plt.grid(linestyle='--')
# plt.ylim(-2, 2)
# plt.xlim(-2, 2)
# plt.show()

plot_output(real_ak, False, 'asda')
plot_output(real_signal_imp, False, 'asda')
plot_output(out_equ_x, False, 'asda')

#
#plot_constellations(noise_on, r1_quad, r2_conj_quad, r1_skew, r2_conj_skew, r1_impairments, r2_conj_impairments):
plot_constellations(noise,' polarization X ', r1_x_quad, r2_x_conj_quad, r1_x_skew, r2_x_conj_skew, r1_x_impairments, r2_x_conj_impairments )
# plot_constellations(noise,' polarization Y ', r1_y_quad, r2_y_conj_quad, r1_y_skew, r2_y_conj_skew, r1_y_impairments, r2_y_conj_impairments )


# # print("Channel 1 symbols for pol X: ", ak_x_ch1)
# # print('')
# # print("Channel 1 symbols for pol Y: ", bk_y_ch1)
# # print('')
# # print("Channel 2 symbols for pol X: ", ak_x_ch2)
# # print('')
# # print("Channel 2 symbols for pol Y: ", bk_y_ch2)
# # print('')
# # print('')
# # print('')
# #
# # print("Rx for channel 1, pol X: ", r1_x_impairments)
# # print('')
# # print("Rx conj for channel 1, pol Y: ", r1_y_impairments)
# # print('')
# # print("Rx for channel 2, conjugated for pol X: ", r2_x_conj_impairments)
# # print('')
# # print("Rx for channel 2, conjugated for pol Y: ", r2_y_conj_impairments)
# # print('')
#
#
# '''    Formato para entradas al MIMO de 8x8:   r_pol_channel_phase/quad(I/Q)
# '''
#
# #Entradas Eq MIMO
# r_X_I_1 = np.real(r1_x_impairments)
# r_X_Q_1 = np.imag(r1_x_impairments)
# r_Y_I_1 = np.real(r1_y_impairments)
# r_Y_Q_1 = np.imag(r1_y_impairments)
#
# r_X_I_2 = np.real(np.conj(r2_x_conj_impairments))
# r_X_Q_2 = np.imag(np.conj(r2_x_conj_impairments))
# r_Y_I_2 = np.real(np.conj(r2_y_conj_impairments))
# r_Y_Q_2 = np.imag(np.conj(r2_y_conj_impairments))
#
#
#
# '''Creo filtros con respuesta al impulso, y con coeficientes iguales a ceros para armar matriz'''
# impulse = np.zeros(filter_taps, dtype = float)
# impulse[int(filter_taps/2)] = 1
#
# zero_filter = np.zeros(int(filter_taps), dtype=float)
#
#
#
#
