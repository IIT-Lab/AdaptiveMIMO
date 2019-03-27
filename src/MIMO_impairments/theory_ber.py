from scipy.special import erfc, binom
import numpy as np

def qam_ber(M, EbN0):
    """ computes the Bit Error Ratio for M-QAM for AWGN as given in
        "Digital Communications", 2nd Ed., Sklar and Ray, [9.54]
    """
    M = np.sqrt(M)
    Pb = (1-M**(-1))/np.log2(M)*erfc(np.sqrt((3*np.log2(M)/(M**2-1)*EbN0)))
    return Pb

def qam_mrc_ber_hsnr(M, L, EbN0):
    """ computes the Bit Error Ratio for M-QAM in the case of Rayleigh
        fading, AWGN and MRC with L branches as given in
        "Digital Communications", 2nd Ed., Sklar and Ray, [9.54]
        "Digital Communications", 4th Ed., Proakis, [14.4-15]
    """
    gamma = EbN0
    d = 1.5*np.log2(M)/(M-1)
    mu = np.sqrt((d*gamma)/(1+d*gamma))
    tot = 0
    for k in np.arange(L):
        tot = tot + binom(L-1+k, k)*np.power((1+mu)/2, k)
    SM = np.sqrt(M)
    P2 = 2*(SM-1)/(SM*np.log2(SM))*np.power((1-mu)/2, L)*tot
    return P2

def qam_mrc_ber_loglin(M, L, EbN0):
    """ computes the Bit Error Ratio for M-QAM in the case of Rayleigh
        fading, AWGN and MRC with L branches as given in
        "Digital Communications", 2nd Ed., Sklar and Ray, [9.54]
        "Digital Communications", 4th Ed., Proakis, [14.4-18]
    """
    gamma = EbN0
    d = 6*np.log2(M)/(M-1)
    SM = np.sqrt(M)
    P2 = 2*(SM-1)/(SM*np.log2(SM))*np.power(1./(d*EbN0), L)*binom(2*L-1, L)
    return P2

def qam_mrc_ber(M, L, EbN0):
    """ computes the Bit Error Ratio for M-QAM in the case of Rayleigh
        fading, AWGN and MRC with L branches as in
        http://www.mathworks.de/de/help/comm/ug/bit-error-rate-ber.html#bq5fugo
        which itself uses:
        "On the General BER Expression of One- and Two-Deimensional Amplitude
        Modulations", Chon and Yoon, IEEE Transactions on Communications, Vol.
        50, 2002.
        "A Unified Approch to the Performance Analysis of Digital Communication
        over Generalized Fading Chanels", Simon and Alouini, Proceedings of the
        IEEE, Vol. 86, 1998
    """
    Pb = 2/(np.pi*np.sqrt(M)*np.log2(np.sqrt(M)))
    tot = 0;
    for k in np.arange(1, np.log2(np.sqrt(M)+1)):
        for i in np.arange((1-2**(-k))*np.sqrt(M)):
            t1 = (-1)**np.floor(i*2**(k-1)/np.sqrt(M))
            t2 = 2**(k-1)-np.floor(i*2**(k-1)/np.sqrt(M)+1/2.)
            # numerical integration for the integral term
            gamma = np.log2(M)*EbN0
            Mg = lambda s: 1/(1-s*gamma)
            thetas = np.linspace(0.001, np.pi/2, 300)
            C = 3*((2*i+1)**2)/(2*(M-1))
            y = np.array([ Mg(-C/((np.sin(theta))**2)) for theta in thetas ])
            y = np.power(y, L)
            t3 = np.trapz(y, thetas)
            tot = tot+t1*t2*t3
    return Pb*tot


# import matplotlib as mp
# mp.rc('text', usetex=True)
# mp.rc('font', family='serif', size=12)
# import matplotlib.pyplot as plt
# # 16-QAM modulation
# M = 16
# # number of receive antennas
# N = np.arange(1, 5)
# # noise vector
# EbN0 = np.logspace(-2, 2, 41)
# # linestyles
# ls = ['b', 'r', 'g', 'k']
# ax = plt.subplot(111)
#
# for idx, n in enumerate(N):
#     bers = np.array([qam_mrc_ber_loglin(M, n, x) for x in EbN0])
#     plt.semilogy(10*np.log10(EbN0), bers, '--' + ls[idx], label='approx. loglin ' + str(n) + ' paths')
#
#     bers = np.array([qam_mrc_ber_hsnr(M, n, x) for x in EbN0])
#     plt.semilogy(10*np.log10(EbN0), bers, '.' + ls[idx], label='approx. high SNR ' + str(n) + ' paths')
#
#     bers = np.array([qam_mrc_ber(M, n, x) for x in EbN0])
#     plt.semilogy(10*np.log10(EbN0), bers, ls[idx], label= 'true ' + str(n) + ' paths')
#
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels, loc=3)
#
# x1,x2,y1,y2 = plt.axis()
# plt.axis((x1,x2,y1,1))
#
# plt.title(str(M) + "-QAM MRC")
# plt.ylabel(r'BER')
# plt.xlabel(r'$E_b/N_0$ in dB')
#
# plt.grid()
# plt.show()