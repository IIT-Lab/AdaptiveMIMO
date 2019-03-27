from configurable import Configurable
import scipy.signal as sign
import numpy as np

class TestChannel(Configurable):
    def __init__(self, ConfigClass):
        self.complex_type = True
        self.taps = 16
        self.brate = 32e9 # [Samples/s]
        self.Fs = self.brate * 1.1
        self.tscale, self.coefficients = TestChannel.sinc(self.taps, self.brate, self.Fs)
        self.coefficients = np.expand_dims(self.coefficients, axis=0)
        self.SNR = None
        self.load_config(ConfigClass)

    def process(self, xn):
        if self.SNR:
            x_power = np.var(xn, axis=1, ddof=1)
            sigma = np.sqrt(10**(-self.SNR/10)*np.sum(x_power)/xn.shape[0])
            noise = sigma*np.random.randn(*xn.shape)
            if self.complex_type:
                noise = noise + 1j*sigma*np.random.randn(*xn.shape)
        else:
            noise = 0
        y = sign.convolve(xn, self.coefficients, mode='same')
        return y + noise

    @staticmethod
    def sinc(N, Brate, Srate, symetric=True):
        t_aux = np.linspace(0.0, N/Srate, N, endpoint=False)
        t = (t_aux - N/(2*Srate)) if symetric else t_aux
        return t, np.sinc(2*Brate*t)
