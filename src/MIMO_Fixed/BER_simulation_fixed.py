import numpy as np
from scipy.special import erfc
#from scipy import weave
#from scipy.weave import converters
import matplotlib.pyplot as plt
from time import time
# from tools._fixedInt import *
from scipy.special import erfc
from DSP_utils import *
from adaptive_LMS_equalizer import LMS_equalizer as LMSequ
from adaptive_LMS_equalizer import LMS_equalizer_fixed as LMS_fixed
from tool import _fixedInt as fix
import funciones






SNRdB       = 12.0 # valores de SNR (en dB) a utilizar
noise_on    = True
channel     = True
n_symbols   = 500000
os_tx       = 1
#fbaud       = 32e9
Tbaud       = 1.0/1024000.0
Nbauds_tx   = 9.0
beta_tx     = 0.51

'''Parametros para representacion en punto fijo '''
NBT = 14
NBF = 12






def get_BER_errors_float(signal, sent_symbols, modulation):

    if (str(modulation))=='QPSK':
        symbols_ref = np.array([[-1. + 1j * 1.], [-1. - 1j * 1.], [1. + 1j * 1.0], [1.0 - 1j * 1.0]])
    elif (str(modulation))=='BPSK':
        symbols_ref = np.array([[-1.0], [1.0]])


    output_slicer = np.zeros(len(signal), dtype='float64')

    for sample in range(0, len(signal)):
        argmin = np.argmin((np.abs(signal[sample] - symbols_ref)) ** 2)

        output_slicer[sample] = symbols_ref[argmin, 0]


    errors = float(np.sum(sent_symbols!=output_slicer))/float(len(sent_symbols))

    return errors

def get_BER_errors_fixed(signal, sent_symbols, modulation):

    if (str(modulation))=='QPSK':
        symbols_ref = np.array([[-1. + 1j * 1.], [-1. - 1j * 1.], [1. + 1j * 1.0], [1.0 - 1j * 1.0]])
    elif (str(modulation))=='BPSK':
        symbols_ref = np.array([[-1.0], [1.0]])


    output_slicer = np.zeros(len(signal), dtype='float64')
    for sample in range(0, len(signal)):

        argmin = np.argmin((np.abs(signal[sample] - symbols_ref)) ** 2)
        output_slicer[sample] = symbols_ref[argmin, 0]


    errors = float(np.sum(sent_symbols!=output_slicer))/float(len(sent_symbols))
    # slicer = funciones.fix_append(fix.arrayFixedInt(NBT,NBF,output_slicer,'S','round','saturate'),0)
    # sent_symbols = funciones.fix_append(fix.arrayFixedInt(NBT,NBF,sent_symbols,'S','round','saturate'),0)

    #errors = float(np.sum(sent_symbols != output_slicer)) / float(len(sent_symbols))

    return errors








#step = 0.0001
equalizer_taps = 21
len_filter = np.copy(equalizer_taps)
#step = 0.0001
step = 0.001
delay_equalizer = int((len_filter-1)/2)
# Creamos instancia del ecualizador
#LMS_equalizer = LMSequ(equalizer_taps, step)
LMS_equalizer_FX = LMS_fixed(equalizer_taps, step)
LMS_equalizer_FLOAT = LMSequ(equalizer_taps, step)

bpsk_symbols = generate_symbols(2,n_symbols)
Ea = np.mean(np.abs(bpsk_symbols)**2)
#
# print(bpsk_symbols)

#Parametros para calcular BER
Eb_N0_dB = np.arange(1,11,1) #-3 a 10 dB





sim_BER_BPSK_float     = []
sim_BER_BPSK_equ_float = []

sim_BER_BPSK_fx       = []
sim_BER_BPSK_equ_fx    = []


terrors          = []

#(t, rrc_tx) = rcosine(beta_tx, Tbaud, os_tx, Nbauds_tx, Norm=True)
#delay_rc = int((len(rrc_tx) - 1) / 2)



for EbN0 in Eb_N0_dB:

    #out_signal = np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,bpsk_symbols,'S','round','saturate'), 0))
    out_signal = np.copy(bpsk_symbols)

    if channel:
        #hk = [0.04, -0.05, 0.07, -0.21, -0.5, 0.72, 0.36, 0, 0.21, 0.03, 0.07]  # anda
        #hk = [0., 1., 0.]
        #hk = [0.407, 0.815, 0.407]
        hk = [-0.0041, 0.0091, -0.0235, 0.0465, -0.0723, 0.0926, 0.9034, 0.0926, -0.0723, 0.0465, -0.0235, 0.0091, -0.0041]
        #hk = [0.05, -0.063, 0.088, -0.126, -0.25, 0.9047, 0.25, 0.0, 0.126, 0.038, 0.088]
        hk = np.array(hk, dtype='float64')
        '''Energia canal'''
        Eh_total = 0
        for var in hk:
            Eh_total = Eh_total + (abs(var)) ** 2
        Eh = np.sqrt(Eh_total)
        hk = hk / float(Eh)

        channel_fixed = np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,hk,'S','round','saturate'), 0))
        #print('La diferencia entre el canal float menos el fix es', np.array(hk - channel_fixed))
        delay_channel = int((len(hk) - 1) / 2)

        out_channel_float = np.convolve(hk, out_signal)
        out_channel_fixed = funciones.fix_append(fix.arrayFixedInt(NBT, NBF,np.convolve(channel_fixed, out_signal),'S','round','saturate'),0)
        out_channel_fixed = reshape_filter_output(out_channel_fixed, delay_channel)
        out_channel_float = reshape_filter_output(out_channel_float, delay_channel)

        #out_channel_fixed = np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,out_channel_fixed,'S','round','saturate'), 0))


        out_float = np.copy(out_channel_float)
        out_fixed =np.copy(out_channel_fixed)


    else:
        out_float = np.copy(out_signal)
        out_fixed= np.copy(out_signal)
        Eh = 1


    if noise_on:
        #awgn = generate_AWGN(out_channel,EbN0,False, Eh, Ea )]

        constellation = np.array([-1., +1.])
        Energy_Bpsk = np.sqrt(np.mean(np.abs(constellation) ** 2))  # Average energy
        EbNo = 10 ** (0.1 * EbN0)
        varianza = Energy_Bpsk / (2.0 * EbNo) # 2*
        DesvEstandar = np.sqrt(varianza)
        awgn = DesvEstandar * np.array(np.random.randn(1,n_symbols)[0])

        rx_signal = out_float+ awgn
        rx_signal_fixed = funciones.fix_append(fix.arrayFixedInt(NBT, NBF,rx_signal,'S','round','saturate'),0)



    else:
        rx_signal = out_float
        rx_signal_fixed = out_fixed

    '''Instancia de ecualizador en punto flotante: recibe la senal representada en punto flotante'''
    out_eq_float, err_eq_float= LMS_equalizer_FLOAT.adapt(rx_signal)

    '''Intsancia de ecualizador en punto flotante: recibe la senal representada en punto fijo'''
    out_eq_fx, err_eq_fx = LMS_equalizer_FX.adapt(rx_signal_fixed)


    #Hacemos el reshape por el largo del ecualizador
    reshape_symbols = bpsk_symbols[delay_equalizer:len(bpsk_symbols)-(delay_equalizer+1)]

    '''Obtencion de errores'''
    error_no_equalized_FLOAT = get_BER_errors_float(rx_signal, bpsk_symbols, 'BPSK')
    error_equalized_FLOAT = get_BER_errors_float(out_eq_float, bpsk_symbols[delay_equalizer:len(bpsk_symbols)-(delay_equalizer+1)], 'BPSK')

    error_no_equalized_FX = get_BER_errors_fixed(rx_signal, bpsk_symbols, 'BPSK')
    #error_equalized_FX = get_BER_errors_fixed(out_eq_fx, bpsk_symbols[delay_equalizer:len(bpsk_symbols)-(delay_equalizer+1)], 'BPSK')
    error_equalized_FX = get_BER_errors_fixed(out_eq_fx, bpsk_symbols[delay_equalizer:len(bpsk_symbols) - (delay_equalizer + 1)], 'BPSK')

    # error_no_equalized_FX = get_BER_errors(rx_signal_fixed, bpsk_symbols, 'BPSK')
    # error_equalized_FX = get_BER_errors(out_eq_fx, bpsk_symbols[delay_equalizer:len(bpsk_symbols) - (delay_equalizer + 1)], 'BPSK')
    sim_BER_BPSK_float.append(error_no_equalized_FLOAT)
    sim_BER_BPSK_equ_float.append(error_equalized_FLOAT)

    sim_BER_BPSK_fx.append(error_no_equalized_FX)
    sim_BER_BPSK_equ_fx.append(error_equalized_FX)

    teorico = float(0.5 * erfc(np.sqrt(10 ** (EbN0 / 10.0))))
    terrors.append(teorico)

    print('SNR', EbN0)
    print('Error teorico: ', teorico)
    print('Error calculado sin EQU FLOAT: ', error_no_equalized_FLOAT)
    print('Error calculado con EQU FLOAT: ', error_equalized_FLOAT)
    print('Error calculado sin EQU FX: ', error_no_equalized_FX)
    print('Error calculado con EQU FX: ', error_equalized_FX)
    print('')

#
print(terrors)
print('')
print('Float')
print(sim_BER_BPSK_float)
print('')
print(sim_BER_BPSK_equ_float)
print('')
print('Fixed')
print(sim_BER_BPSK_fx)
print('')
print(sim_BER_BPSK_equ_fx)

# #
# sim_BER_BPSK_fx = np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,sim_BER_BPSK_fx,'S','round','saturate'), 0))
# sim_BER_BPSK_equ_fx= np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,sim_BER_BPSK_equ_fx,'S','round','saturate'), 0))
# #
# #
#





#
f = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=10)
plt.semilogy(Eb_N0_dB,sim_BER_BPSK_float,'-o', label = 'No Equalized Signal(ISI + Noise) FLOAT')
plt.semilogy(Eb_N0_dB,sim_BER_BPSK_equ_float,'-x', label = 'Equalized Signal FLOAT')
#plt.semilogy(Eb_N0_dB,sim_BER_BPSK_fx,'-o', label = 'No Equalized Signal(ISI + Noise) FIXED')
plt.semilogy(Eb_N0_dB,sim_BER_BPSK_equ_fx,'-s', label = 'Equalized Signal FIXED')
#plt.semilogy(Eb_N0_dB,sim_BER_Q_NOISE_noEq,'-x', label = 'No equalized: Signal + Noise (No channel response)')
plt.semilogy(Eb_N0_dB,terrors,':v', label='Theoritical')
#plt.legend(('simulation','theory'),'upper right',shadow=False,fancybox=True)
plt.title('Bit Error Rate')
plt.xlabel(r'$\frac{Eb}{N0}$, dB')
plt.ylabel('Error probability')
plt.legend()
plt.grid(True)
plt.savefig('Plots/' + str('BER_curve_fixed_channel3_more_samples') + '.png')
plt.show()




file= open("BER_results.txt","a")
file.write('Resultado de simulacion para NBT = ' +str(NBT) + '  y NBF = ' + str(NBF) + 'con N_symbols = ' + str(n_symbols))
file.write('\n')
file.write('Vector BER teorico:  ' + str(terrors))
file.write('\n')
file.write('Vector BER para senal no ecualizada en punto flotante:   ' + str(sim_BER_BPSK_float))
file.write('\n')
file.write('Vector BER para senal ecualizada en punto flotante:   ' + str(sim_BER_BPSK_equ_float))
file.write('\n')
file.write('Vector BER para senal no ecualizada en punto fijo:   ' + str(sim_BER_BPSK_fx))
file.write('\n')
file.write('Vector BER para senal ecualizada en punto fijo:   ' + str(sim_BER_BPSK_equ_fx))
file.write('\n\n')
file.close()

