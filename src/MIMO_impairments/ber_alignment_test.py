import numpy as np
from scipy.special import erfc
#from scipy import weave
#from scipy.weave import converters
import matplotlib.pyplot as plt
from time import time
# from tools._fixedInt import *
from tools.DSPtools import*
from tools.r_cosine import*
from scipy.special import erfc
from DSP_utils import *
from adaptive_LMS_equalizer import LMS_equalizer as LMSequ


SNRdB       = 12.0 # valores de SNR (en dB) a utilizar
noise_on    = False
channel     = True
RC          = False
n_symbols   = 100
os_tx       = 1
#fbaud       = 32e9
Tbaud       = 1.0/1024000.0
Nbauds_tx   = 11.0
beta_tx     = 0.51
N = 12000


def slicer_and_error(complex_signal, sent_I, sent_Q, N):


 ipHat_I = complex_signal.real>0
 ipHat_Q = complex_signal.imag>0
 err_I = float((sent_I != ipHat_I).sum())/N
 err_Q = float((sent_Q != ipHat_Q).sum()) / N

 return err_I, err_Q










#step = 0.0001
equalizer_taps = 15
len_filter = np.copy(equalizer_taps)
step = 0.004

delay_equalizer = int((len_filter-1)/2)



#Creamos instancia del ecualizador
LMS_equalizer = LMSequ(equalizer_taps, step)

#Generamos simbolos
ip_I = np.random.rand(1,N)>0.5;
ch1_I = 2*ip_I-1
ip_Q = np.random.rand(1,N)>0.5;
ch1_Q = 2*ip_Q-1
complex_ch1 = ch1_I + (1j*ch1_Q)
n = 1.0/np.sqrt(2.0)*(np.random.randn(1,N)+1j*np.random.randn(1,N))



#Parametros para calcular BER
Eb_N0_dB = np.arange(1.0,14.0,1.0) #-3 a 10 dB



sim_BER_I_NOISE_noEq= []
sim_BER_Q_NOISE_noEq = []

sim_BER_I_NOISE_Eq = []
sim_BER_Q_NOISE_Eq = []



terrors = []




for EbN0 in Eb_N0_dB:
 if RC:
  #filtro caida cosenoidal (no srrc), i.e, simulamos Tx y Rx en cascada (normalizado)
  (t, rrc_tx) = rcosine(beta_tx, Tbaud, os_tx, Nbauds_tx, Norm=True)
  delay_rc = int((len(rrc_tx)-1)/2)
  out_signal = np.convolve(rrc_tx, complex_ch1[0])
  out_signal = reshape_filter_output(out_signal, delay_rc)
  out_signal = np.reshape(out_signal, (1, len(out_signal)))
 else:
  out_signal = complex_ch1

 if channel:
  #hk = [0.04, -0.05, 0.07, -0.21, -0.5, 0.72, 0.36, 0, 0.21, 0.03, 0.07]  # anda
  # hk con 13 coef
  hk = [-0.0041, 0.0091, -0.0235, 0.0465, -0.0723, 0.0926, 0.9034, 0.0926, -0.0723, 0.0465, -0.0235, 0.0091, -0.0041]
  hk = np.array(hk)
  '''Energia canal'''
  Eh = np.sum(hk ** 2)
  hk = np.array(hk) / np.sqrt(Eh)

  out_channel = np.convolve(hk, out_signal[0])
  delay_channel = int((len(hk) - 1) / 2)
  out_channel = reshape_filter_output(out_channel, delay_channel)

  out_signal = np.reshape(out_channel, (1, len(out_channel)))


 rx_signal = np.copy(out_signal)

 if noise_on:
  rx_signal = out_signal + (10 ** (-EbN0 / 20.0) * n)




 sent_signal_I = np.copy(rx_signal.real)
 sent_signal_Q = np.copy(rx_signal.imag)


 output_eq_I, error_out_I = LMS_equalizer.adapt(sent_signal_I[0])
 output_eq_Q, error_out_Q = LMS_equalizer.adapt(sent_signal_Q[0])


 out_eq_complex = output_eq_I + (1j*output_eq_Q)
 out_eq_complex = np.reshape(out_eq_complex, (1,len(out_eq_complex)))

 sent_I = ip_I[0,delay_equalizer:(len(ip_I[0])-(delay_equalizer+1))]
 sent_Q = ip_Q[0, delay_equalizer:(len(ip_Q[0])-(delay_equalizer+1))]
 #print(ip_Q)


 '''Error simulado con unicamente AWGN'''
 err_eq_I, err_eq_Q = slicer_and_error(out_eq_complex, sent_I, sent_Q, N)
 err_I_noise, err_Q_noise = slicer_and_error(rx_signal, ip_I, ip_Q, N)


 sim_BER_I_NOISE_noEq.append(err_I_noise)
 sim_BER_Q_NOISE_noEq.append(err_Q_noise)
 sim_BER_I_NOISE_Eq.append(err_eq_I)
 sim_BER_Q_NOISE_Eq.append(err_eq_Q)

 terrors.append(float(0.5*erfc(np.sqrt(10**(EbN0/10.0)))))
 #print("time = %s, EbN0 = %s, ber = %s"%(time()-t,EbN0,err_I))




print(terrors)
print('')
print(sim_BER_I_NOISE_noEq)
print(sim_BER_Q_NOISE_noEq)
print('')
print(sim_BER_I_NOISE_Eq)
print(sim_BER_Q_NOISE_Eq)

# #f = plt.figure()
# file_name = 'BER_curves_I'
# plt.figure(figsize=(10, 8))
# plt.title('Bit Error Rate for LMS Equalizer with step ' + str(step))
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=10)
# plt.semilogy(Eb_N0_dB,sim_BER_I_NOISE_noEq,'-o', label = 'No Equalized Signal(Signal + Noise)')
# plt.semilogy(Eb_N0_dB,sim_BER_I_NOISE_Eq,'-x', label = 'Equalized Signal (Signal + Noise)')
# #plt.semilogy(Eb_N0_dB,sim_BER_Q_NOISE_noEq,'-x', label = 'No equalized: Signal + Noise (No channel response)')
# plt.semilogy(Eb_N0_dB,terrors,':v', label='Theoritical')
# #plt.legend(('simulation','theory'),'upper right',shadow=False,fancybox=True)
# plt.title('Bit Error Rate')
# plt.xlabel(r'$\frac{Eb}{N0}$, dB')
# plt.ylabel('Bit Error Rate')
# plt.legend()
# plt.grid(True)
# plt.savefig('Plots/' + str(file_name) + '.pdf')
# plt.show()
#
#
#

f = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=10)
plt.semilogy(Eb_N0_dB,sim_BER_I_NOISE_noEq,'-o', label = 'No Equalized Signal(Signal + Noise)')
plt.semilogy(Eb_N0_dB,sim_BER_I_NOISE_Eq,'-x', label = 'Equalized Signal (Signal + Noise)')
#plt.semilogy(Eb_N0_dB,sim_BER_Q_NOISE_noEq,'-x', label = 'No equalized: Signal + Noise (No channel response)')
plt.semilogy(Eb_N0_dB,terrors,':v', label='Theoritical')
#plt.legend(('simulation','theory'),'upper right',shadow=False,fancybox=True)
plt.title('Bit Error Rate')
plt.xlabel(r'$\frac{Eb}{N0}$, dB')
plt.ylabel('Bit Error Rate')
plt.legend()
plt.grid(True)
plt.show()


plot_output(output_eq_I,False, ';asjfkasdjsda')



