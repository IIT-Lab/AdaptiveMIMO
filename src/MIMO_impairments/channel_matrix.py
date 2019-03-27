#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy             as np
import matplotlib.pyplot as plt
import scipy.signal      as signal
import scipy.signal      as linalg
import scipy.special     as special
import math
import sys
import os

from cmath import *
from decimal import Decimal    as ToDec



#########################
# Funciones y Variables #
#########################
# Number of Bits to transmit
n_data = 5000
# Bandwidth [Baud]
fBaud = 32e9
# Baud rate
Tbaud = 1 / fBaud
# Upsampling factor
OS_factor = 8
# Sample rate
Ts = Tbaud / OS_factor
# Bandwidth [Hz]
Fs = 1 / Ts
# Roll Off factor
roll_off = 0.50
# Offset bandwith
freq_offset = (1 + roll_off) / (2 * Tbaud)
# Turn noise on or off
switch_channel = True * 1

# Range for BER curve
SNR_range = np.arange(1, 10)
Perror = np.zeros(len(SNR_range), dtype='float')
# Ber counters for analysis
ber_I = np.arange(len(SNR_range), dtype='float') * 0
ber_Q = np.arange(len(SNR_range), dtype='float') * 0
# Ber counters CASO 3
ber_I_1xt = np.arange(len(SNR_range), dtype='float') * 0
ber_Q_2xc = np.arange(len(SNR_range), dtype='float') * 0

####################
# Simbolos Canal 1 #
####################
simbolos_ch1 = (2 * np.random.randint(2, size=(2, n_data))) - 1
ch1_I = simbolos_ch1[0,]
ch1_Q = simbolos_ch1[1,]
ch1 = ch1_I + 1j * ch1_Q

####################
# Simbolos Canal 2 #
####################
simbolos_ch2 = (2 * np.random.randint(2, size=(2, n_data))) - 1
ch2_I = simbolos_ch2[0,]
ch2_Q = simbolos_ch2[1,]
ch2 = ch2_I + 1j * ch2_Q

###########################
# CASO 1 Error Quadrature #
###########################
# tau=0T ,eg=0.2, phi=20°
eg = 0.1 * 1
phi_e = 20 * 1
phi_e = (phi_e * np.pi) / 180.0

Alfa = 0.5 * ((1 - eg) * np.exp(1j * phi_e / 2) + (1 + eg) * np.exp(-1j * phi_e / 2))
Beta = 0.5 * ((1 - eg) * np.exp(-1j * phi_e / 2) - (1 + eg) * np.exp(+1j * phi_e / 2))
M_phi_eg = np.array([[Alfa, Beta], [np.conj(Beta), np.conj(Alfa)]])

# Matriz de Simbolos complejos
a1n = (ch1)
a2n = (ch2)
a12 = np.array([a1n, np.conj(a2n)])

# Receptor con Efectos fase y amplitud
r12 = np.matmul(M_phi_eg, a12)
r1x = r12[0,]
r2xc = r12[1,]

# # Visualizacion
# print('Parte real de r1x'	, np.real(r1x)[:15])
# print( 'Parte Imaginaria de r1x'	, np.imag(r1x)[:15])
# print('Parte real de r2xc'		, np.real(r2xc)[:15])
# print('Parte Imaginaria de r2xc', np.imag(r2xc)[:15])
# Constellation
plt.figure ();
plt.title  ('Receptor: Symbols + Error Quadrature');
# plt.scatter(np.real(Alfa*a1n + Beta*np.conj(a1n)), np.imag(Alfa*a1n + Beta*np.conj(a1n)), marker='*', c='r');#np.imag(r1x)
plt.scatter(np.real(r1x), np.imag(r2xc), marker='.', c='lime')  ;  # np.imag(r1x)
plt.grid(linestyle='--');
plt.ylim(-2, 2); plt.xlim(-2, 2)


###############
# CASO 2 Skew #
###############
# tau=0.1T ,eg=0, phi=0
# Si |tau|<0.2T, la Matriz es valida
# Es satisfactorio pára Qpsk, sin embargo, se degrada para Qam16
tau 	 = 0.0625 *Tbaud *1  # 0.0625
wc  	 = 2* np.pi * freq_offset
M_tau = np.array(
    [[np.cos(wc * tau / 2), 1j * np.sin(wc * tau / 2)], [+1j * np.sin(wc * tau / 2), np.cos(wc * tau / 2)]])
# Receptor con Efectos Skew
r12_tau = np.matmul(M_tau, a12)
r1x_tau = r12_tau[0,]
r2xc_tau = r12_tau[1,]

# Visualizacion
# print'\nParte real de r1x_tau'	, np.real(r1x_tau)[:15]
# print 'Parte Imaginaria de r1x_tau'	, np.imag(r1x_tau)[:15]

plt.figure ()
plt.title  ('Receptor: Symbols + Skew')
plt.scatter(np.real(r1x_tau), np.imag(r2xc_tau), marker='.', c='fuchsia')
plt.grid(linestyle='--')
plt.ylim(-2, 2); plt.xlim(-2, 2)


############################
# CASO 3 Skew + Quadrature #
############################
M_t    = np.matmul(M_tau, M_phi_eg)
# Receptor Total
r12_t  = np.matmul(M_t, a12)
r1x_t  = r12_t[0,]
r2xc_t = r12_t[1,]

# # Visualizacion
# print '\nParte real de r1x_t'		, np.real(r1x_t)[:15]
# print 'Parte Imaginaria de r1x_t'	, np.imag(r1x_t)[:15]
plt.figure ();
plt.title  ('Receptor: Efecto Completo');
plt.scatter(np.real(r1x_t),np.ones(len(r1x_t)) ,marker='.', c='r');
plt.grid(linestyle='--');
plt.ylim(-2, 2); plt.xlim(-2, 2)


################
# BER Simulado #
################
for ptr, dB_values in enumerate(SNR_range):
    # for ptr, dB_values in enumerate([10]):
    ##############
    # Ruido AWGN #
    ##############
    if switch_channel:
        SNR_dB   	  = dB_values  # 8
        constellation = np.array([ -1. +1.j, -1. -1.j, 1. +1.j,  1. -1.j])
        constellation = np.array([ -1., +1.])
        Energy_Bpsk   = np.sqrt(np.mean(np.abs(constellation )**2)  )  # Average energy
        EbNo	 	  = 10 **(0.1 *SNR_dB)
        varianza 	  = Energy_Bpsk /(2.0 * EbNo)  # 2*
        awgn_I    	  = np.sqrt(varianza) * np.random.randn(1, n_data)[0]
        awgn_Q 	 	  = np.sqrt(varianza) * np.random.randn(1, n_data)[0]

        awgn      	  = awgn_I + 1j *awgn_Q
        ch1_noise 	  = ch1 +awgn
        ch2_noise 	  = ch2 +awgn

        # Caso 3
        r1x_t_noise_ber  = r1x_t  + awgn
        r2xc_t_noise_ber = r2xc_t + awgn

    else:
        awgn 	  = ( 0 +0j ) *np.arange(n_data)
        ch1_noise = ch1 +awgn
        ch2_noise = ch2 +awgn

        # Caso 3
        r1x_t_noise_ber  = r1x_t  + awgn
        r2xc_t_noise_ber = r2xc_t + awgn


    ############
    # Decision #
    ############
    # Caso 3
    r1x_t_hat 	  = [0 ] *len(r1x_t)
    r2xc_t_hat 	  = [0 ] *len(r2xc_t)
    # Canal sin efectos
    ch1_noise_hat = [0 ] *len(ch1)
    symbols_ref   = np.array([ [-1. + 1j *1.], [-1. - 1j *1.], [1. + 1j *1.0], [1.0 - 1j *1.0] ])

    for muestra in range(len(ch1)):
        argmin        		   = np.argmin(  (np.abs(ch1_noise[muestra] - symbols_ref) )**2  )
        ch1_noise_hat[muestra] = symbols_ref[argmin ,0]

        # Caso 3
        argmin_1x        	   = np.argmin(  (np.abs(r1x_t_noise_ber[muestra] - symbols_ref) )**2  )
        r1x_t_hat[muestra] 	   = symbols_ref[argmin_1x ,0]
        argmin_2xc        	   = np.argmin(  (np.abs(r2xc_t_noise_ber[muestra] - symbols_ref) )**2  )
        r2xc_t_hat[muestra]    = symbols_ref[argmin_2xc ,0]


    ##########################
    # Contador Errores Canal #
    ##########################
    for i ,q in enumerate(range(n_data)):
        if np.real(ch1_noise_hat)[i] != np.real(ch1)[i] :
            ber_I[ptr] = ber_I[ptr] + 1

        if np.imag(ch1_noise_hat)[q] != np.imag(ch1)[q] :
            ber_Q[ptr] = ber_Q[ptr] + 1

        # Caso 3
        if np.real(r1x_t_hat)[i] != np.real(a1n)[i] :
            ber_I_1xt[ptr] = ber_I_1xt[ptr] + 1

        if np.imag(r2xc_t_hat)[q] != np.imag(np.conj(a2n))[q] :
            ber_Q_2xc[ptr] = ber_Q_2xc[ptr] + 1

    # print 'SNRdB=%4.2f, BER I=%10.4e, BER Q=%10.4e' % (SNR_range[ptr], ber_I[ptr],ber_Q[ptr])
    #print 'SNRdB=%4.2f, BER I=%10.4e, BER Q=%10.4e' % (SNR_range[ptr], ber_I_1xt[ptr], ber_Q_2xc[ptr])


##################
# BER Simulacion #
##################
Ber_iq = []
Ber_12 = []
for value in range(len(SNR_range)):
    Ber_iq.append( (ber_I[value] + ber_Q[value]) / (1.0 * n_data) )
    # Caso 3
    Ber_12.append( (ber_I_1xt[value] + ber_Q_2xc[value]) / (1.0 * n_data) )

##########################
# BER Teorico 			 #
##########################
# Para calcular Ber agregar Ruido y Alinear las secuencias Tx*Rx
# Gaussian Q-function
def Q_function(x):
    return 0.5 *special.erfc( x /np.sqrt(2.))

ProbabilityE_s = np.arange(len(SNR_range), dtype = 'float' ) *0
ProbabilityE_b = np.arange(len(SNR_range), dtype = 'float' ) *0

for ptr, dB_values in enumerate(SNR_range):
    SNR_dB 		        = dB_values
    EbNo	 	        = 10**(0.1 *SNR_dB)
    ProbabilityE_s[ptr] = 2* Q_function(np.sqrt(2 * EbNo)) * (1 - 0.5 * Q_function(np.sqrt(2 * EbNo)))
    ProbabilityE_b[ptr] = Q_function(np.sqrt(2 * EbNo))
# print 'PE_theoretical = %10.4e' % ProbabilityE_s[ptr]

# Visualizacion
plt.figure();
plt.title('QPSK - Curva de SER')
plt.semilogy(SNR_range, ProbabilityE_s, 'r-', label="Curva Teorica")
plt.semilogy(SNR_range, Ber_iq, 'm+:', label="Curva simulacion: Canal+Ruido")
plt.semilogy(SNR_range, Ber_12, 'g*:', label="Curva simulacion: Canal+efectos+Ruido")

plt.xlabel('Signal to Noise Ratio [dB]');
plt.ylabel('Error Probability');
plt.legend(loc=0);
plt.minorticks_on();
plt.grid(which='minor', linestyle='--')

########################
# Constellation Caso 1 #
########################
r1x_noise = r1x + awgn
r2xc_noise = r2xc + awgn

r1x_tau_noise = r1x_tau + awgn
r2xc_tau_noise = r2xc_tau + awgn

r1x_t_noise = r1x_t + awgn
r2xc_t_noise = r2xc_t + awgn




plt.figure();
plt.title('QPSK Constelacion: Error AmplitudFase + Ruido')
plt.scatter(np.real(r1x_noise), np.imag(r2xc_noise), marker='.', c='lime');
plt.grid(linestyle='--');
plt.ylim(-2, 2);
plt.xlim(-2, 2)

########################
# Constellation Caso 2 #
########################

plt.figure();
plt.title('QPSK Constelacion: Error Skew + Ruido')
plt.scatter(np.real(r1x_tau_noise), np.imag(r2xc_tau_noise), marker='.', c='fuchsia');
plt.grid(linestyle='--');
plt.ylim(-2, 2);
plt.xlim(-2, 2)

########################
# Constellation Caso 3 #
########################

plt.figure();
plt.title('QPSK Constelacion: Error AmplitudFase + Skew + Ruido')
plt.scatter(np.real(r1x_t_noise), np.imag(r2xc_t_noise), marker='.', c='r');
plt.grid(linestyle='--');
plt.ylim(-2, 2);
plt.xlim(-2, 2)


plt.figure();
plt.title('QPSK Constelacion canal 1')
plt.scatter(np.real(r1x_t_noise), np.imag(r1x_t_noise), marker='.', c='r');
plt.grid(linestyle='--');
plt.ylim(-2, 2);
plt.xlim(-2, 2)


plt.figure();
plt.title('QPSK Constelacion canal 2')
plt.scatter(np.real(r2xc_t_noise), np.imag(r2xc_t_noise), marker='.', c='r');
plt.grid(linestyle='--');
plt.ylim(-2, 2);
plt.xlim(-2, 2)


plt.figure();
plt.title('QPSK Constelacion al reves')
plt.scatter(np.real(r2xc_t_noise), np.imag(r1x_t_noise), marker='.', c='r');
plt.grid(linestyle='--');
plt.ylim(-2, 2);
plt.xlim(-2, 2)

###################################
plt.show()
########## ##########################




print('')
print(len(ch1))
print(' ')
print(len(ch2))
print(' ')
print(len(r1x))
print(' ')
print(len(r2xc))
print(' ')
print(len(r1x_tau))
print(' ')
print(len(r2xc_tau))
print(' ')
print(len(r1x_t))
print(' ')
print(len(r2xc_t))
print(' ')
print(len(r1x_noise))
print(' ')
print(len(r2xc_noise))
print(' ')
print(len(r1x_tau_noise))
print(' ')
print(len(r2xc_tau_noise))
print(' ')
print(len(r1x_t_noise))
print(' ')
print(len(r2xc_t_noise))
print(' ')




















