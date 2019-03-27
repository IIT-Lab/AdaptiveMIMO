import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
# from tools._fixedInt import *
from tools.DSPtools import*
from tools.r_cosine import*
from scipy.special import erfc
from DSP_utils import *

from adaptive_LMS_equalizer import LMS_equalizer as LMSequ


SNRdB       = 12.0 # valores de SNR (en dB) a utilizar
noise_on    = True
channel     = True
SRRC        = False
n_symbols   = 100000
os_tx       = 1
fbaud       = 32e9
Tbaud       = 1/fbaud
Nbauds_tx   = 8.0
beta_tx     = 0.5



# # Sample rate
# Ts = Tbaud / os_tx
# # Bandwidth [Hz]
# Fs = 1 / Ts


#step = 0.0001
equalizer_taps = 21
len_filter = np.copy(equalizer_taps)
step = 0.0004
delay = int((len_filter-1)/2)






#Creamos instancia del ecualizador
LMS_equalizer = LMSequ(equalizer_taps, step)

'''Generacion de simbolos'''
qpsk_symbols, symbols_I, symbols_Q = generate_symbols(M = 4, n_symbols = n_symbols)


#qpsk_symbols = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
#
# qpsk_symbols = np.ones(n_symbols)
# qpsk_symbols[0] = 0
# qpsk_symbols[len(qpsk_symbols)-1] = 0
symbols_I = np.ones(n_symbols)
symbols_I[0] = 0
symbols_I[len(symbols_I)-1] = 0
symbols_Q = np.ones(n_symbols)
symbols_Q[0] = 0
symbols_Q[len(symbols_Q)-1] = 0


Ea = np.mean(np.abs(symbols_I)**2)
print(Ea)

'''Conformacion de pulso con SRRC'''
if SRRC:
    (t, rrc_tx) = rcosine(beta_tx, Tbaud, os_tx, Nbauds_tx, Norm=False)  # Filtro Tx en punto flotante
    out_tx = np.convolve(rrc_tx, qpsk_symbols)
    # same convolution with even length arrays introduces a 0 in 1st pos
    out_tx = out_tx[int((Nbauds_tx - 1) / 2) + 1:int(len(out_tx)) - int((Nbauds_tx - 1) / 2)]
else:
    out_tx = np.copy(qpsk_symbols)


if channel:
    #hk = [0.04, -0.05, 0.07, -0.21, -0.5, 0.72, 0.36, 0, 0.21, 0.03, 0.07]  # anda
    #hk = [-0.0041, 0.0091, -0.0235, 0.0465, -0.0723, 0.0926, 0.9034, 0.0926, -0.0723, 0.0465, -0.0235, 0.0091,-0.0041]
    hk = [0.05, -0.063, 0.088, -0.126, -0.25, 0.9047, 0.25, 0.0, 0.126, 0.038, 0.088]
    #hk = [-0.0041, 0.0091, -0.0235, 0.0465, -0.0723, 0.0926, 0.9034, 0.0926, -0.0723, 0.0465, -0.0235, 0.0091, -0.0041]
    hk = np.array(hk)
    '''Energia canal'''
    # Eh_total = 0
    # for var in hk:
    #     Eh_total = Eh_total + (abs(var)) ** 2
    #     Eh = np.sqrt(Eh_total)
    #Normalizo canal
    Eh = np.sum(hk ** 2)
    hk = np.array(hk)/np.sqrt(Eh)

    #Convolucion senal transmitida con el canal
    out_channel = np.convolve(out_tx, hk)

    #Reshape por la convolucion que agrega elementos
    out_channel = out_channel[int(len(hk) / 2) - 1:len(out_channel) - int((len(hk) - 1) / 2)]



    # La convolucion con el impulso, no agrega nada.

else:
    out_channel = np.copy(out_tx)
    print(out_channel)
    Eh = 1

'''Genero AWGN'''
if noise_on ==True:

    noise = generate_AWGN(out_channel,SNRdB, True, Eh, Ea)
    sent_signal = np.copy(out_channel) + noise


    # noise = generate_AWGN(out_channel, SNRdB, complex)
    # print('Ruido es:', noise)
    # print('Mande hasta ahora:', out_channel)
    # sent_signal = np.copy(out_channel) + noise
    # print('Suma es:', sent_signal)

else:
    sent_signal = np.copy(out_channel)

#Separo lineas I/Q
sent_signal_I = np.real(sent_signal)
#sent_signal_Q = np.imag(sent_signal)


#steps = [0.115, 0.045, 0.009, 0.0001, 0.00001, 0.000001, 0.025, 0.004]

output_eq_I, error_out_I = LMS_equalizer.adapt(sent_signal_I)
#output_eq_Q, error_out_Q = LMS_equalizer.adapt(sent_signal_Q)

'''Slicer'''
output_slicer_I= np.sign(output_eq_I)
#output_slicer_Q= np.sign(output_eq_Q)

#plot_constellation(sent_signal_I, sent_signal_Q, False, 'Constelacion de salida')
# #plot_constellation(output_eq_I, output_eq_Q, False, 'Constelacion de salida')
#
# plt.figure()
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=10)
# plt.title('Constelacion de salida')
# #plt.plot(sent_signal_I,np.zeros(len(sent_signal_I)),'.', linewidth = 2.0)
# plt.plot(output_eq_I,np.zeros(len(output_eq_I)),'.', linewidth = 2.0)
# #plt.plot(sent_signal_I, sent_signal_Q, 'x', linewidth=2.0)
# #plt.plot(output_eq_I, output_eq_Q, '.', linewidth = 2.0)
# plt.xlabel('Parte real')
# plt.ylabel('Parte imaginaria')
# plt.grid(True)
# plt.show()



# for i in range(0, 100):
#     plt.figure()
#     plt.stem(np.array(LMS_equalizer.c_history[i]))
#     plt.show()


#
#
# #Debo detectar desde symbols_I[delay:len(symbols_I)-delay-1]
# #Esto se debe a que al pasar por un filtro, los n taps del mismo te dan el delay
# #En cuanto al filtro del canal, y el SRRC, ya hice el reshape antes, por eso no lo tengo en cuenta ahora,
# #si en el ecualizador
#
# #print(np.array_equal(output_slicer, symbols_I[delay:len(symbols_I)-delay-1]))
#
# # #print(symbols_I[delay:len(symbols_I)-delay-1])
# #
# plot_impulse_response(LMS_equalizer.coefficients,False)
# #plot_impulse_response(LMS_equalizer.coefficients, True)
# plot_output(output_eq_I, False,'Output Equalizer for step = ' + str(step))
# plot_output(sent_signal_I, False,'Senal enviada')
#plot_output(output_eq_I, True, 'Output Equalizer for step = ' + str(step))
# signals_to_plot = [hk, LMS_equalizer.coefficients, np.convolve(hk,LMS_equalizer.coefficients)]
# plot_freq_responses(signals_to_plot, ['Channel Response', 'Equalizer Response', 'Output'], True,'Frequency responses')
















#print(len(LMS_equalizer.c_history))
#plot_output(output_eq_I, True, 'Output Equalizer for step ' + str(LMS_equalizer.step))

# fig = plt.figure(figsize=(10, 8))
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=10)
# ax1 = plt.subplot(221)
# plt.title("iteracion 1000")
# plt.stem(LMS_equalizer.c_history[1000])
# plt.grid()
# ax2 = plt.subplot(222, sharey=ax1)
# plt.title("h12")
# plt.title("iteracion 2000")
# plt.stem(LMS_equalizer.c_history[2000])
# plt.grid()
# ax3 = plt.subplot(223, sharex=ax1)
# plt.title("iteracion 100000")
# plt.stem(LMS_equalizer.c_history[10000])
# plt.grid()
# ax4 = plt.subplot(224, sharey=ax3, sharex=ax1)
# plt.title("Ultima iteracion ")
# plt.stem(LMS_equalizer.c_history[99900])
# plt.grid()
# plt.savefig('Plots/' + str('Coeficientes') + '.png')
# plt.show()

#
# print(LMS_equalizer.c_history)
#
# plt.figure()
# plt.title('Channel response')
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=10)
# plt.stem(hk)
# plt.grid()
# plt.xlabel('n')
# plt.ylabel('Amplitude')
# plt.savefig('Plots/' + str('out_channel_coef2') + '.png')
# plt.show()

#
# plot_output(output_eq_I, True, 'Output equalizer for step ' + str(LMS_equalizer.step))
# plot_output(sent_signal_I, True, 'Signal without equalization')