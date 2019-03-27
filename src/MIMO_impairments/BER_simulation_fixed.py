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
import commpy.filters as compyf
from tools import _fixedInt as fix






SNRdB       = 12.0 # valores de SNR (en dB) a utilizar
noise_on    = True
channel     = True
RC          = False
n_symbols   = 5000
os_tx       = 1
#fbaud       = 32e9
Tbaud       = 1.0/1024000.0
Nbauds_tx   = 9.0
beta_tx     = 0.51

'''Parametros para representacion en punto fijo '''
NBT = 14
NBF = 12






def get_BER_errors(signal, sent_symbols, modulation):

    if (str(modulation))=='QPSK':
        symbols_ref = np.array([[-1. + 1j * 1.], [-1. - 1j * 1.], [1. + 1j * 1.0], [1.0 - 1j * 1.0]])
    elif (str(modulation))=='BPSK':
        symbols_ref = np.array([[-1.0], [1.0]])

    signal = fix.arrayFixedInt(NBT,NBF,signal,'S','round','saturate')
    sent_symbols = fix.arrayFixedInt(NBT,NBF,sent_symbols,'S','round','saturate')
    output_slicer = np.zeros(len(signal), dtype='float64')

    for sample in range(0, len(signal)):
        # print('')
        # print('Recibi', signal[sample])
        # print('Mande', sent_symbols[sample])
        argmin = np.argmin((np.abs(signal[sample] - symbols_ref)) ** 2)
        #argmin = np.argmin((np.abs(signal[sample] - sent_symbols[sample])) ** 2)
        # print('Detecte', symbols_ref[argmin, 0])
        output_slicer[sample] = symbols_ref[argmin, 0]


    errors = float(np.sum(sent_symbols!=output_slicer))/float(len(sent_symbols))
    errors = fix.arrayFixedInt(NBT,NBF,errors,'S','round','saturate')
    #print('Error calculado:', str(errors))
    return errors







#step = 0.0001
equalizer_taps = 21
len_filter = np.copy(equalizer_taps)
#step = 0.0001
step = 0.0001
delay_equalizer = int((len_filter-1)/2)
# Creamos instancia del ecualizador
#LMS_equalizer = LMSequ(equalizer_taps, step)
LMS_equalizer = LMS_fixed(equalizer_taps, step)


bpsk_symbols = generate_symbols(2,n_symbols)
Ea = np.mean(np.abs(bpsk_symbols)**2)
#
# print(bpsk_symbols)

#Parametros para calcular BER
Eb_N0_dB = np.arange(1,10,1) #-3 a 10 dB





sim_BER_BPSK     = []
sim_BER_BPSK_equ = []
terrors          = []

#(t, rrc_tx) = rcosine(beta_tx, Tbaud, os_tx, Nbauds_tx, Norm=True)
#delay_rc = int((len(rrc_tx) - 1) / 2)



for EbN0 in Eb_N0_dB:

    out_signal = fix.arrayFixedInt(NBT,NBF,bpsk_symbols,'S','round','saturate')

    if channel:
        #hk = [0.04, -0.05, 0.07, -0.21, -0.5, 0.72, 0.36, 0, 0.21, 0.03, 0.07]  # anda
        #hk = [0., 1., 0.]
        #hk = [0.407, 0.815, 0.407]
        #hk = [-0.0041, 0.0091, -0.0235, 0.0465, -0.0723, 0.0926, 0.9034, 0.0926, -0.0723, 0.0465, -0.0235, 0.0091, -0.0041]
        hk = [0.05, -0.063, 0.088, -0.126, -0.25, 0.9047, 0.25, 0.0, 0.126, 0.038, 0.088]
        hk = np.array(hk, dtype='float64')
        '''Energia canal'''
        Eh_total = 0
        for var in hk:
            Eh_total = Eh_total + (abs(var)) ** 2
        Eh = np.sqrt(Eh_total)
        hk = hk / float(Eh)

        out_channel = np.convolve(hk, out_signal)
        delay_channel = int((len(hk) - 1) / 2)
        out_channel = reshape_filter_output(out_channel, delay_channel)

    else:
        out_channel = np.copy(out_signal)
        Eh = 1


    if noise_on:
        #awgn = generate_AWGN(out_channel,EbN0,False, Eh, Ea )]

        constellation = np.array([-1., +1.])
        Energy_Bpsk = np.sqrt(np.mean(np.abs(constellation) ** 2))  # Average energy
        EbNo = float(10 ** (0.1 * EbN0))
        varianza = float(Energy_Bpsk / (2.0 * EbNo))  # 2*
        DesvEstandar = float(np.sqrt(varianza))
        awgn = DesvEstandar * np.array(np.random.randn(1,n_symbols)[0], dtype='float64')
        rx_signal = out_channel + awgn
    else:
        rx_signal = np.copy(out_channel)


    #step = float(1/(5*equalizer_taps*np.var(rx_signal,ddof=1)))
    #LMS_equalizer.step = step

    out_eq, err_eq= LMS_equalizer.adapt(rx_signal)


    #Perdes los primeros y los ultimos
    reshape_symbols = bpsk_symbols[delay_equalizer:len(bpsk_symbols)-(delay_equalizer+1)]

    # print(len(detected_equ))
    # print(len(detected_no_equ))


    # out_slicer_equ = np.sign(out_eq)
    # out_slicer_NoEqu = np.sign(rx_signal)
    # #

    error_no_equalized = get_BER_errors(rx_signal, bpsk_symbols, 'BPSK')
    error_equalized = get_BER_errors(out_eq, bpsk_symbols[delay_equalizer:len(bpsk_symbols)-(delay_equalizer+1)], 'BPSK')



    sim_BER_BPSK.append(error_no_equalized)
    sim_BER_BPSK_equ.append(error_equalized)
    teorico = float(0.5 * erfc(np.sqrt(10 ** (EbN0 / 10.0))))
    terrors.append(teorico)

    print('SNR', EbN0)
    print('Error teorico: ', teorico)
    print('Error calculado sin EQU: ', error_no_equalized )
    print('Error calculado con EQU: ', error_equalized)
    print('')


print(terrors)
print('')
print(sim_BER_BPSK)
print('')
print(sim_BER_BPSK_equ)

#
f = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=10)
plt.semilogy(Eb_N0_dB,sim_BER_BPSK,'-o', label = 'No Equalized Signal(ISI + Noise)')
plt.semilogy(Eb_N0_dB,sim_BER_BPSK_equ,'-x', label = 'Equalized Signal')
#plt.semilogy(Eb_N0_dB,sim_BER_Q_NOISE_noEq,'-x', label = 'No equalized: Signal + Noise (No channel response)')
plt.semilogy(Eb_N0_dB,terrors,':v', label='Theoritical')
#plt.legend(('simulation','theory'),'upper right',shadow=False,fancybox=True)
plt.title('Bit Error Rate')
plt.xlabel(r'$\frac{Eb}{N0}$, dB')
plt.ylabel('Error probability')
plt.legend()
plt.grid(True)
#plt.savefig('Plots/' + str('BER_curve') + '.png')
plt.show()

# mse_obtained = np.array(LMS_equalizer.mse)
# mse_obtained = mse_obtained/len(out_eq)
# print(mse_obtained)
# plt.figure()
# #plt.plot(mse_eq)
# plt.semilogy(mse_obtained)
# #plt.plot(mse_obtained)
# #plt.semilogy(np.mean(err_eq))
# #plt.semilogy(np.abs(err_eq))
# #plt.plt(10*np.log)
# plt.show()

