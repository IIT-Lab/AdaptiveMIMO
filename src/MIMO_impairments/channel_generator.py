#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy             as np
import matplotlib.pyplot as plt
import scipy.signal      as signal
import scipy.signal      as linalg
import scipy.special     as special
import math
import sys

from cmath import *
from decimal import Decimal    as ToDec
from tools._fixedInt import *
from tools import DSPtools   as dsp



'''Variables para almacenar errores'''
# Range for BER curve
SNR_range = np.arange(1, 11)  # (1,10)#(10,19)
Perror = np.zeros(len(SNR_range), dtype='float')
# Ber counters for analysis
ber_I = np.arange(len(SNR_range), dtype='float') * 0
ber_Q = np.arange(len(SNR_range), dtype='float') * 0
# Ber counters CASO 3
ber_I_1xt = np.arange(len(SNR_range), dtype='float') * 0
ber_Q_2xc = np.arange(len(SNR_range), dtype='float') * 0
# Caso CUANTIZACION
ber_I_1xt_c = np.arange(len(SNR_range), dtype='float') * 0
ber_Q_2xc_c = np.arange(len(SNR_range), dtype='float') * 0



'''Parametros del sistema'''

dual_polarization = False
# Number of Bits to transmit
n_data = 80000
# Number of carriers
n_scarriers = 2
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


'''Parametros con efectos a incluir'''
#Skew
tau = 0.05 * Tbaud * 1  # 0.0625
#Error ganancia
eg = 0.2 * 1
#Error de fase
phi_e = 12 * 1
phi_e = (phi_e * np.pi) / 180.0
#Frecuencia de portadora (para modulacion)
wc = 2 * np.pi * freq_offset


'''Parametros para representacion en punto fijo '''
NBT = 14
NBF = 12












'''Generacion simbolos canal 1'''
#Genero matriz de 2x2 con simbolos aleatorios +-1
symbols_ch1 = (-2 * np.random.randint(2, size=(2, n_data))) + 1
#Primer vector asignado a componente en fase
ch1_I = symbols_ch1[0,]
#Segundo vector o fila asignado a componente en cuadratura
ch1_Q = symbols_ch1[1,]
#Vector con simbolos complejos correspondiente al primer canal
ch1_complex = ch1_I + 1j * ch1_Q


'''Generacion simbolos canal 2'''
#Genero matriz de 2x2 con simbolos aleatorios +-1
symbols_ch2 = (-2 * np.random.randint(2, size=(2, n_data))) + 1
#Primer vector asignado a componente en fase
ch2_I = symbols_ch2[0,]
#Segundo vector o fila asignado a componente en cuadratura
ch2_Q = symbols_ch2[1,]
#Vector con simbolos complejos correspondiente al primer canal
ch2_complex = ch2_I + 1j * ch2_Q



'''Componentes de los simbolos'''
a1I = ch1_I
a1Q = ch1_Q

a2I = ch2_I
a2Q = ch2_Q


'''Constantes a CUANTIZAR segun la formula, incluyendo los impairments'''
A = np.cos(wc * tau / 2 + phi_e / 2)
B = np.cos(wc * tau / 2 - phi_e / 2)
C = np.sin(wc * tau / 2 + phi_e / 2)
D = np.sin(wc * tau / 2 - phi_e / 2)
E = eg * np.cos(wc * tau / 2 + phi_e / 2)
F = eg * np.cos(wc * tau / 2 - phi_e / 2)
G = eg * np.sin(wc * tau / 2 + phi_e / 2)
H = eg * np.sin(wc * tau / 2 - phi_e / 2)


'''Componentes de las señales recibidas con impairments'''
#Senal recibida primer canal (r1x)
Real1x = a1I * B + a1Q * G + a2I * (-E) + a2Q * D
Imag1x = a1I * (-G) + a1Q * B + a2I * (D) + a2Q * E
#Senal recibida segundo canal CONJUGADA (r2* = r2xc)
Real2xc = a1I * (-F) + a1Q * (-C) + a2I * A + a2Q * (-H)
Imag2xc = a1I * C + a1Q * (-F) + a2I * (-H) + a2Q * (-A)

r1x_complex = Real1x + 1j * Imag1x
r2xc_complex = Real2xc + 1j * Imag2xc

# # Visualizacion Parametros
# print
# 'Parametros Ideales:\n'
# print
# 'A = ', A
# print
# 'B = ', B
# print
# 'C = ', C
# print
# 'D = ', D
# print
# 'E = ', E
# print
# 'F = ', F
# print
# 'G = ', G
# print
# 'H = ', H
# print
# '\n'


#Constellation Caso 3 |Sin| Ruido
# plt.figure(1)
# plt.title('Received signal with impairments included (without noise)');
# plt.scatter(np.real(r1x_complex), np.imag(r1x_complex), marker='.', c='lime')  # np.imag(r1x)
# plt.grid(linestyle='--')
# plt.ylim(-2, 2)
# plt.xlim(-2, 2)
# plt.show()
#
# plt.figure(2)
# plt.title('Received signal with impairments included (without noise)');
# plt.scatter(np.real(r2xc_complex), np.imag(r2xc_complex), marker='.', c='lime')  # np.imag(r1x)
# plt.grid(linestyle='--')
# plt.ylim(-2, 2)
# plt.xlim(-2, 2)
# plt.show()
#
#
# plt.figure(3)
# plt.title('Received signal with impairments included (without noise)');
# plt.scatter(np.real(r1x_complex), np.imag(r2xc_complex), marker='.', c='lime')  # np.imag(r1x)
# plt.grid(linestyle='--')
# plt.ylim(-2, 2)
# plt.xlim(-2, 2)
# plt.show()
#
# plt.figure(4)
# plt.title('Received signal with impairments included (without noise)');
# plt.scatter(np.real(r2xc_complex), np.imag(r1x_complex), marker='.', c='lime')  # np.imag(r1x)
# plt.grid(linestyle='--')
# plt.ylim(-2, 2)
# plt.xlim(-2, 2)
# plt.show()






# ################
#
# Parametros = [A, B, C, D, E, F, G, H]
# Parametros_fp = arrayFixedInt(NBT, NBF, Parametros, 'S', 'round', 'saturate')  # 16,14
# Parametros_C = [0] * len(Parametros_fp)
# for ptr in range(len(Parametros_fp)):
#     Parametros_C[ptr] = Parametros_fp[ptr].fValue
#     print
#     'Parametros Cuantizados:', Parametros_C[ptr], '| Int:%d |' % Parametros_fp[ptr].intvalue, 'Bin:', bin(
#         Parametros_fp[ptr].intvalue)
#
# # Componentes CUANTIZADAS
# A_c = Parametros_C[0]
# B_c = Parametros_C[1]
# C_c = Parametros_C[2]
# D_c = Parametros_C[3]
# E_c = Parametros_C[4]
# F_c = Parametros_C[5]
# G_c = Parametros_C[6]
# H_c = Parametros_C[7]
#
# Real1x_c = a1I * B_c + a1Q * G_c + a2I * (-E_c) + a2Q * D_c
# Imag1x_c = a1I * (-G_c) + a1Q * B_c + a2I * (D_c) + a2Q * E_c
#
# Real2xc_c = a1I * (-F_c) + a1Q * (-C_c) + a2I * A_c + a2Q * (-H_c)
# Imag2xc_c = a1I * C_c + a1Q * (-F_c) + a2I * (-H_c) + a2Q * (-A_c)
#
# R1x_c = Real1x_c + 1j * Imag1x_c
# R2xc_c = Real2xc_c + 1j * Imag2xc_c
#
# # for ptr in range(20):
# # 	print 'Real1x_c: ',Real1x_c[ptr], '\tch1_I: ',ch1_I[ptr]
#
################
# BER Simulado #
################
for ptr, dB_values in enumerate(SNR_range):
    # for ptr, dB_values in enumerate([15]):
    ##############
    # Ruido AWGN #
    ##############
    if switch_channel:
        SNR_dB = dB_values  # 8
        constellation = np.array([-1. + 1.j, -1. - 1.j, 1. + 1.j, 1. - 1.j])
        constellation = np.array([-1., +1.])
        Energy_Bpsk = np.sqrt(np.mean(np.abs(constellation) ** 2))  # Average energy
        EbNo = 10 ** (0.1 * SNR_dB)
        varianza = Energy_Bpsk / (2.0 * EbNo)  # 2*
        DesvEstandar = np.sqrt(varianza)
        awgn_I = DesvEstandar * np.random.randn(1, n_data)[0]
        awgn_Q = DesvEstandar * np.random.randn(1, n_data)[0]

        awgn = awgn_I + 1j * awgn_Q
        ch1_noise = ch1_complex + awgn
        ch2_noise = ch2_complex + awgn

        # Caso 3 IDEAL
        r1x_t_noise_ber = r1x_complex + awgn
        r2xc_t_noise_ber = r2xc_complex + awgn

#         ################################
#         # CUANTIZACION Desvio Estandar #
#         ################################
#         DE_bits = 10
#         DE_frac = 8
#         DesvEstandar_c = DeFixedInt(DE_bits, DE_frac, 'S', 'round', 'saturate')
#         DesvEstandar_c.value = DesvEstandar
#         print
#         'SNRdB:', SNR_range[
#             ptr], '| DE:', DesvEstandar, '| DE_c:', DesvEstandar_c.fValue, '| Int:%d|' % DesvEstandar_c.intvalue, 'Bin:', bin(
#             DesvEstandar_c.intvalue)
#         print
#         '\n'
#
#         # CUANTIZACION del Ruido
#         R_bit = 16
#         R_frac = 14
#         Ruido_I_fp = arrayFixedInt(R_bit, R_frac, awgn_I, 'S', 'round', 'saturate')
#         Ruido_Q_fp = arrayFixedInt(R_bit, R_frac, awgn_Q, 'S', 'round', 'saturate')
#
#         awgn_I_c = [0] * len(Ruido_I_fp)
#         awgn_Q_c = [0] * len(Ruido_I_fp)
#
#         for n in range(len(Ruido_I_fp)):
#             awgn_I_c[n] = Ruido_I_fp[n].fValue
#             awgn_Q_c[n] = Ruido_Q_fp[n].fValue
#
#         awgn_c = np.array(awgn_I_c) + 1j * np.array(awgn_Q_c)
#
#         # Caso señal CUANTIZADA
#         R1x_c_noise = R1x_c + awgn_c
#         R2xc_c_noise = R2xc_c + awgn_c
#
    else:
        awgn = (0 + 0j) * np.arange(n_data)
        ch1_noise = ch1_complex + awgn
        ch2_noise = ch2_complex + awgn

        # Caso 3 IDEAL
        r1x_t_noise_ber = r1x_complex + awgn
        r2xc_t_noise_ber = r2xc_complex + awgn
#         # Caso CUANTIZADO
#         R1x_c_noise = R1x_c + awgn
#         R2xc_c_noise = R2xc_c + awgn
#
#     ################################################
#     # Cuantizo la salida Total Señal+Efectos+Ruido #
#     ################################################
#     RT_bit = 16  # Range: -4.0000000000 up to 3.9998779297
#     RT_frac = 13
#
#     RT1_I_fp = arrayFixedInt(RT_bit, RT_frac, np.real(R1x_c_noise), 'S', 'round', 'saturate')
#     RT1_Q_fp = arrayFixedInt(RT_bit, RT_frac, np.imag(R1x_c_noise), 'S', 'round', 'saturate')
#     RT2_I_fp = arrayFixedInt(RT_bit, RT_frac, np.real(R2xc_c_noise), 'S', 'round', 'saturate')
#     RT2_Q_fp = arrayFixedInt(RT_bit, RT_frac, np.imag(R2xc_c_noise), 'S', 'round', 'saturate')
#
#     RT1_I_c = [0] * len(RT1_I_fp)
#     RT1_Q_c = [0] * len(RT1_Q_fp)
#     RT2_I_c = [0] * len(RT2_I_fp)
#     RT2_Q_c = [0] * len(RT2_Q_fp)
#
#     for n in range(len(RT1_I_fp)):
#         RT1_I_c[n] = RT1_I_fp[n].fValue
#         RT1_Q_c[n] = RT1_Q_fp[n].fValue
#         RT2_I_c[n] = RT2_I_fp[n].fValue
#         RT2_Q_c[n] = RT2_Q_fp[n].fValue
#
#     R1x_c_noise = np.array(RT1_I_c) + 1j * np.array(RT1_Q_c)
#     R2xc_c_noise = np.array(RT2_I_c) + 1j * np.array(RT2_Q_c)
#
#     ############
#     # Decision #
#     ############
    # Caso 3 IDEAL
    r1x_t_hat = [0] * len(r1x_complex)
    r2xc_t_hat = [0] * len(r2xc_complex)
#     # Caso CUANTIZADO
#     R1x_c_hat = [0] * len(R1x_c)
#     R2xc_c_hat = [0] * len(R2xc_c)
#
    # Canal sin efectos
    ch1_noise_hat = [0] * len(ch1_complex)
    symbols_ref = np.array([[-1. + 1j * 1.], [-1. - 1j * 1.], [1. + 1j * 1.0], [1.0 - 1j * 1.0]])

    for muestra in range(len(ch1_complex)):
        argmin = np.argmin((np.abs(ch1_noise[muestra] - symbols_ref)) ** 2)
        ch1_noise_hat[muestra] = symbols_ref[argmin, 0]

        # Caso 3 IDEAL
        argmin_1x = np.argmin((np.abs(r1x_t_noise_ber[muestra] - symbols_ref)) ** 2)
        r1x_t_hat[muestra] = symbols_ref[argmin_1x, 0]
        argmin_2xc = np.argmin((np.abs(r2xc_t_noise_ber[muestra] - symbols_ref)) ** 2)
        r2xc_t_hat[muestra] = symbols_ref[argmin_2xc, 0]
#
#         # Caso CUANTIZADO
#         argmin_1x_c = np.argmin((np.abs(R1x_c_noise[muestra] - symbols_ref)) ** 2)
#         R1x_c_hat[muestra] = symbols_ref[argmin_1x_c, 0]
#         argmin_2xc_c = np.argmin((np.abs(R2xc_c_noise[muestra] - symbols_ref)) ** 2)
#         R2xc_c_hat[muestra] = symbols_ref[argmin_2xc_c, 0]
#
    ##########################
    # Contador Errores Canal #
    ##########################
    real1x = np.real(r1x_t_hat)
    imag2xc = np.imag(r2xc_t_hat)
#     real1x_c = np.real(R1x_c_hat)
#     imag2xc_c = np.imag(R2xc_c_hat)
#
    real_ch1 = np.real(ch1_complex)
    imagC_ch2 = np.imag(np.conj(ch2_complex))
#
    for i in range(n_data):
        q = i
        # Caso 3 IDEAL
        if real1x[i] != real_ch1[i]:
            ber_I_1xt[ptr] = ber_I_1xt[ptr] + 1

        if imag2xc[q] != imagC_ch2[q]:
            ber_Q_2xc[ptr] = ber_Q_2xc[ptr] + 1
#
#         # Caso CUANTIZADO
#         if real1x_c[i] != real_ch1[i]:
#             ber_I_1xt_c[ptr] = ber_I_1xt_c[ptr] + 1
#
#         if imag2xc_c[q] != imagC_ch2[q]:
#             ber_Q_2xc_c[ptr] = ber_Q_2xc_c[ptr] + 1
#
#     # print 'SNRdB=%4.2f, BER I=%10.4e, BER Q=%10.4e' % (SNR_range[ptr], ber_I[ptr],ber_Q[ptr])
#     print
#     'SNRdB=%4.2f, BER I=%10.4e, BER Q=%10.4e' % (SNR_range[ptr], ber_I_1xt_c[ptr], ber_Q_2xc_c[ptr])
#
##################
# BER Simulacion #
##################
Ber_12 = []
# Ber_c = []
for value in range(len(SNR_range)):
    # Caso 3 IDEAL
    Ber_12.append((ber_I_1xt[value] + ber_Q_2xc[value]) / (2.0 * n_data))
#     # Caso CUANTIZADO
#     Ber_c.append((ber_I_1xt_c[value] + ber_Q_2xc_c[value]) / (2.0 * n_data))
#
#
# ##########################
# # BER Teorico 			 #
# ##########################
# # Para calcular Ber agregar Ruido y Alinear las secuencias Tx,Rx
#
#

ProbabilityE_s = np.arange(len(SNR_range), dtype='float') * 0
ProbabilityE_b = np.arange(len(SNR_range), dtype='float') * 0
#
for ptr, dB_values in enumerate(SNR_range):
    SNR_dB = dB_values
    EbNo = 10 ** (0.1 * SNR_dB)
    ProbabilityE_s[ptr] = 2 * dsp.Q_function(np.sqrt(2 * EbNo)) * (1 - 0.5 * dsp.Q_function(np.sqrt(2 * EbNo)))
    ProbabilityE_b[ptr] = dsp.Q_function(np.sqrt(2 * EbNo))
# print 'PE_theoretical = %10.4e' % ProbabilityE_s[ptr]


# Visualizacion SIMULACION TEORICO
#f = plt.figure();
plt.figure()
plt.title('BER curve')
plt.semilogy(SNR_range, ProbabilityE_b, 'r-', label="Curva Teorica")
plt.semilogy(SNR_range, Ber_12, 'g*:', label="Curva simulacion: Symbol+Noise+Effects(FullResolution)")
# plt.semilogy(SNR_range, Ber_c, 'm+:', label="Curva simulacion: Symbol+Noise+Effects S(14,12)")
plt.xlabel('Signal to Noise Ratio [dB]');
plt.ylabel('Error Probability');
plt.legend(loc=0);
plt.minorticks_on();
plt.grid(which='minor', linestyle='--');
plt.grid(which='mayor', linestyle='-')
#f.savefig("foo.pdf", bbox_inches='tight')








####################################
# Constellation Caso 3 Ideal Ruido #
####################################
r1x_t_noise = r1x_complex + awgn
r2xc_t_noise = r2xc_complex + awgn
plt.figure();
plt.title('QPSK Constelacion: Error AmplitudFase + Skew + Ruido')
#plt.scatter(np.real(R1x_c_noise), np.imag(R2xc_c_noise), marker='.', c='r');
plt.scatter(np.real(r1x_t_noise), np.imag(r2xc_t_noise), marker='.', c='r');
plt.grid(linestyle='--');
plt.ylim(-2, 2);
plt.xlim(-2, 2)

###################################
plt.show()
print
# '\n\n'
# ###################################