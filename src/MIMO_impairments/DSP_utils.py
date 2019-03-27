import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy.special import erfc



def generate_symbols(M, n_symbols):

    if M ==2:
        # Generamos simbolos +-1
        bpsk_symbols = np.random.choice([-1,+1], n_symbols)
        return bpsk_symbols

    elif M==4:
        # Generamos simbolos complejos
        qpsk_symbols = np.random.choice([-1 + 1j, -1 - 1j, 1 + 1j, 1 - 1j], n_symbols)
        # Separamos lineas
        symbols_I = np.real(qpsk_symbols)
        symbols_Q = np.imag(qpsk_symbols)
        return qpsk_symbols, symbols_I, symbols_Q

def generate_AWGN(signal, SNRdB, complex, Eh, Ea):


    # The variance is the average of the squared deviations
    # from the mean, i.e., var = mean(abs(x - x.mean())**2).

    #Ea = np.mean(np.abs(signal)**2)
    N0 = Ea/SNRdB
    if complex:
        noise = (1/np.sqrt(2))* (np.random.randn(len(signal)) + (1j*np.random.randn(len(signal))))
        noise = float((np.sqrt(Eh*N0) )) *noise
    else:
        noise = (1 / np.sqrt(2)) * (np.array(np.random.randn(len(signal)), dtype='float64'))
        noise = float(( np.sqrt(Eh*N0) )) * noise

    return noise






    # # signal_power = np.var(signal)
    # # print(signal_power)
    # # sigma = np.sqrt(10 ** (-SNRdB / 10) * np.sum(signal_power) / signal.shape[0])
    # # noise = sigma * np.random.randn(signal.shape[0])
    # #
    # # if complex == True:
    # #     noise = noise + 1j * sigma * np.random.randn(signal.shape[0])
    # #
    #
    # # return noise
    #
    # Energy_Bpsk = np.sqrt(np.mean(np.abs(signal) ** 2))  # Average energy
    # EbNo = 10 ** (0.1 * SNRdB)
    # varianza = Energy_Bpsk / (2.0 * EbNo)  # 2*
    # DesvEstandar = np.sqrt(varianza)
    # noise = DesvEstandar * np.random.randn(1, len(signal))[0]
    #
    # if complex:
    #     awgn_I = DesvEstandar * np.random.randn(1, len(signal))[0]
    #     awgn_Q = DesvEstandar * np.random.randn(1, len(signal))[0]
    #     noise = awgn_I + (1j*awgn_Q)
    #
    # return noise

def reshape_filter_output(output, filter_delay):
    reshaped = np.array(np.copy(output[(filter_delay):len(output)-filter_delay]))
    return reshaped



'''Plot functions'''
def plot_impulse_response(signal, save):

    if save:
        file_name = input('Insert file name to save plot: ')
        plt.figure(figsize=(10,8))
        plt.title('Impulse response')
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=10)
        plt.stem(signal)
        plt.grid(True)
        plt.xlabel('n')
        plt.ylabel('Amplitude')
        plt.savefig('Plots/' + str(file_name) +'.pdf')
        plt.show()

    else:
        plt.figure(figsize=(10, 8))
        plt.stem(signal)
        plt.grid(True)
        plt.xlabel('n')
        plt.ylabel('Amplitude')
        plt.show()

def plot_output(signal, save, plot_title):
    title = str(plot_title)

    if save:
        file_name = input('Insert file name to save plot: ')
        plt.figure(figsize=(10,8))
        plt.title(plot_title)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=10)
        plt.plot(signal, '.', linewidth=2.0)
        plt.grid(True)
        plt.xlabel('Samples')
        plt.ylabel('Amplitude')
        plt.savefig('Plots/' + str(file_name) + '.png')
        plt.show()

    else:
        plt.figure(figsize=(10, 8))
        plt.title(plot_title)
        plt.plot(signal, '.', linewidth=2.0)
        plt.grid(True)
        plt.xlabel('Samples')
        plt.ylabel('Amplitude')
        plt.show()

def plot_freq_responses(signal_vector,labels, save, plot_title):

    title = str(plot_title)
    if save:
        file_name = input('Insert file name to save: ')
        plt.figure()
        plt.title(title)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=10)
        it = 0
        for signal in signal_vector:
            print(it)
            print(signal)
            m, l = scipy.signal.freqz(signal)
            plt.plot(m, 20 * np.log10(abs(l)), label=str(labels[it]))
            if it == len(signal_vector):
                break
            else:
                it = it + 1
        plt.legend(prop={'size': 8})
        plt.ylabel('Amplitude [dB]', color='b')
        plt.xlabel('Frequency [rad/sample]')
        plt.grid()
        plt.savefig('Plots/' + str(file_name) + '.png')
        plt.show()



    else:
        plt.figure()
        plt.title(title)
        it = 0
        for signal in signal_vector:
            print(it)
            print(signal)
            m,l= scipy.signal.freqz(signal)
            plt.plot(m, 20 * np.log10(abs(l)), label=str(labels[it]))
            if it==len(signal_vector):
                break
            else:
                it = it+1
        plt.legend(prop={'size': 8})
        plt.ylabel('Amplitude [dB]', color='b')
        plt.xlabel('Frequency [rad/sample]')
        plt.show()

def plot_constellation(signal1, signal2, save, plot_title):
    title = str(plot_title)
    if save:
        file_name = input('Insert file name to save plot: ')
        plt.figure(figsize=(10, 8))
        plt.title(plot_title)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=10)
        plt.plot(signal1,signal2,'.',linewidth=2.0)
        # ##plt.xlim(1000,1500)
        plt.xlabel('Parte real')
        plt.ylabel('Parte imaginaria')
        plt.grid(True)
        plt.savefig('Plots/' + str(file_name) + '.png')
        plt.show()

        # # plt.figure(figsize=(10,8))
        # # plt.plot(yk_I_vec[len(yk_I_vec)-500:], yk_Q_vec[len(yk_Q_vec)-500:],'.',linewidth=2.0)
        # # #plt.xlim(1000,1500)
        # # plt.grid(True)
        # # plt.title('Constelacion de salida')
        # # plt.xlabel('Parte real')
        # # plt.ylabel('Parte imaginaria')
        # # plt.show()

    else:

        plt.figure(figsize=(10, 8))
        plt.title(plot_title)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=10)
        #plt.scatter(signal1, signal2, linewidths=2.0)
        plt.plot(signal1[:1000], signal2[:1000], 'x', linewidth=2.0)
        #plt.plot(signal2, ':v', linewidth=2.0)
        # ##plt.xlim(1000,1500)
        plt.xlabel('Parte real')
        plt.ylabel('Parte imaginaria')
        plt.grid(True)
        plt.show()














