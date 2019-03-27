
import numpy             as np
import matplotlib as plt



class Polarization():
    def __init__(self):

        # Number of Bits to transmit
        self.n_data = 200000
        # Bandwidth [Baud]
        self.fBaud = 32e9
        # Baud rate
        self. Tbaud = 1 / self.fBaud
        # Roll Off factor
        self.roll_off = 0.50
        # Offset bandwith
        self.freq_offset = (1 + self.roll_off) / (2 * self.Tbaud)
        # Turn noise on or off
        self.switch_channel = True * 1

        #De momento trabajamos al baudio
        # Upsampling factor: de momento al baudio
        #self.OS_factor = 8
        # Sample rate
        #self.Ts = self.Tbaud / self.OS_factor
        # Bandwidth [Hz]
        #self.Fs = 1 / self.Ts


        #Tx impairments
        self.eg = 0.1 * 1  #Gain error
        self.phi_e = 20 * 1 #Phase error in grades
        self.phi_e = (self.phi_e * np.pi) / 180.0 #Phase error in radians

        self.Alfa = 0.5 * ((1 - self.eg) * np.exp(1j * self.phi_e / 2) + (1 + self.eg) * np.exp(-1j * self.phi_e / 2))
        self. Beta = 0.5 * ((1 -self. eg) * np.exp(-1j * self.phi_e / 2) - (1 + self.eg) * np.exp(+1j * self.phi_e / 2))

        self.tau = 0.0625 * self.Tbaud * 1  # Skew
        self.wc = 2 * np.pi * self.freq_offset



        #Representacion matricial con los impairments incluidos

        #IQ imbalance matrix
        self.M_phi_eg = np.array(
            [        [self.Alfa, self.Beta], [np.conj(self.Beta), np.conj(self.Alfa)]

            ]

        )

        #Skew matrix
        self.M_tau = np.array(
            [
                [np.cos(self.wc * self.tau / 2), 1j * np.sin(self.wc * self.tau / 2)],
                [+1j * np.sin(self.wc * self.tau / 2), np.cos(self.wc * self.tau / 2)]

            ]

        )

        #Total matrix
        self.impairments_matrix = np.matmul(self.M_tau, self.M_phi_eg)




        '''
        Informacion del canal: para no confundir, los atributos solo constan de r1 y r2 para referenciar que puede
        ser polarizacion V o H (o X/Y) junto con _conj si se trata del conjugado de la senal compleja
        
        '''

        # Simbolos complejos genrados originales en cada canal
        self.symb_ch1_complex = [0] * self.n_data
        self.symb_ch2_complex = [0] * self.n_data

        #Simbolos complejos genrados en cada canal para polarizacion X/Y con ruido
        self.ch1_complex = [0 ] * self.n_data
        self.ch2_complex = [0 ] * self.n_data

        #Senales complejas con error de cuadratura incluido unicamente.
        self.r1_quad = [0 ] * self.n_data
        self.r2_conj_quad = [0 ] * self.n_data

        #Senales complejas con error de skew unicamente incluido
        self.r1_skew = [0 ] * self.n_data
        self.r2_conj_skew = [0 ] * self.n_data

        #Senales complejas con todos los impairments incluidos
        self.r1_impairments = [0 ] * self.n_data
        self.r2_conj_impairments = [0 ] * self.n_data


    def generate_noise(self, SNR_dB):

        constellation = np.array([-1. + 1.j, -1. - 1.j, 1. + 1.j, 1. - 1.j])
        #constellation = np.array([-1., +1.])
        Energy_Bpsk = np.sqrt(np.mean(np.abs(constellation) ** 2))  # Average energy
        EbNo = 10 ** (0.1 * SNR_dB)
        varianza = Energy_Bpsk / (2.0 * EbNo)  # 2*
        awgn_I = np.sqrt(varianza) * np.random.randn(1, self.n_data)[0]
        awgn_Q = np.sqrt(varianza) * np.random.randn(1, self.n_data)[0]
        awgn = awgn_I + 1j * awgn_Q

        return awgn

    def generate_lanes(self, noise_on, SNR_DB):

        simbols_ch1 = (2 * np.random.randint(2, size=(2, self.n_data))) - 1
        ch1_I = simbols_ch1[0,]
        ch1_Q = simbols_ch1[1,]
        ch1_complex = ch1_I + 1j * ch1_Q

        simbols_ch2 = (2 * np.random.randint(2, size=(2, self.n_data))) - 1
        ch2_I = simbols_ch2[0,]
        ch2_Q = simbols_ch2[1,]
        ch2_complex = ch2_I + 1j * ch2_Q

        #Guardo simbolos originales generados
        self.symb_ch1_complex = np.copy(ch1_complex)
        self.symb_ch2_complex = np.copy(ch2_complex)


        # Matriz de Simbolos complejos
        a1n = ch1_complex
        a2n = ch2_complex
        symbols_matrix = np.array([a1n, np.conj(a2n)])

        # Receptor con Efectos fase y amplitud
        rx_IQ_imbalance = np.matmul(self.M_phi_eg,  symbols_matrix)
        r1x_quad = rx_IQ_imbalance[0,]
        r2xc_quad = rx_IQ_imbalance[1,]


        '''Receptor unicamente con efectos de Skew'''
        rx_signal_with_skew = np.matmul(self.M_tau, symbols_matrix)
        r1x_skew = rx_signal_with_skew[0,]
        r2xc_skew = rx_signal_with_skew[1,]

        '''Receptor con todos los efectos incluidos (skew, gain, phase error)'''
        rx_with_impairments = np.matmul(self.impairments_matrix, symbols_matrix)
        r1x_impairments = rx_with_impairments[0,]
        r2xc_impairments = rx_with_impairments[1,]

        if noise_on:
            noise_complex = self.generate_noise(SNR_DB)

            ch1_complex = ch1_complex + noise_complex
            ch2_complex = ch2_complex + noise_complex

            r1x_quad = r1x_quad + noise_complex
            r2xc_quad = r2xc_quad + noise_complex

            r1x_skew = r1x_skew + noise_complex
            r2xc_skew = r2xc_skew + noise_complex

            r1x_impairments = r1x_impairments + noise_complex
            r2xc_impairments = r2xc_impairments + noise_complex

        #Actualizo estado objeto: para no confundir,
        self.ch1_complex = ch1_complex
        self.ch2_complex = ch2_complex
        self.r1_quad = r1x_quad
        self.r2_conj_quad = r2xc_quad
        self.r1_skew = r1x_skew
        self.r2_conj_skew = r2xc_skew
        self.r1_impairments = r1x_impairments
        self.r2_conj_impairments = r2xc_impairments











































