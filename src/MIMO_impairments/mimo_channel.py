from configurable import Configurable
import scipy.signal as sign
import numpy as np

class MIMOChannel(Configurable):
    def __init__(self, ConfigClass):
        self.taps = 16
        self.coefficients = np.zeros((2,2,self.taps))
        self.tscale = np.zeros(self.taps)
        self.SNR = None
        self.os_factor = 1
        self.load_config(ConfigClass)

    def process(self, xn):
        if self.SNR:
            x_power = np.var(xn, axis=1, ddof=1)
            sigma = np.sqrt(10**(-self.SNR/10)*np.sum(x_power)/xn.shape[0])
            noise = sigma*np.random.randn(*xn.shape)
            if xn.dtype == np.complex_: # TODO: ver caso de otro complex type
                noise = noise + 1j*sigma*np.random.randn(*xn.shape)
        else:
            noise = 0

        prefix = int(np.ceil((self.taps-1)/2))
        sufix =  int(np.floor((self.taps-1)/2))
        hpad_data = np.zeros((xn.shape[0], xn.shape[1]+prefix+sufix), dtype=xn.dtype)
        hpad_data[:, prefix:-sufix] = xn
        y = np.array([
            np.squeeze(sign.convolve2d(hpad_data, self.coefficients[0,::-1], mode='valid'), axis=0),
            np.squeeze(sign.convolve2d(hpad_data, self.coefficients[1,::-1], mode='valid'), axis=0)
        ])

        return y + noise

if __name__ == '__main__':
    from symbol_generator import SymbolGenerator
    from test_channel import TestChannel
    from plotable import fft_plot
    import matplotlib.pyplot as plt

    class SGenConfig:
        n_tx_symbs = 1e4
        M = 2
        complex_type = True
        dualpol = True
        os_factor = 1

    class ChnConfig:
        taps = 16
        brate = 32e9
        _, h11 = TestChannel.sinc(taps, brate, 1.1*brate)
        _, h12 = TestChannel.sinc(taps, brate, 1.2*brate)
        _, h21 = TestChannel.sinc(taps, brate, 1.3*brate)
        _, h22 = TestChannel.sinc(taps, brate, 1.4*brate)

        _, jh11 = TestChannel.sinc(taps, brate, 1.15*brate)
        _, jh12 = TestChannel.sinc(taps, brate, 1.25*brate)
        _, jh21 = TestChannel.sinc(taps, brate, 1.35*brate)
        _, jh22 = TestChannel.sinc(taps, brate, 1.45*brate)
        coefficients = np.array([
            [h11+1j*jh11, h12+1j*h12],
            [h21+1j*jh21, h22+1j*h22]
        ])

    gen = SymbolGenerator(SGenConfig)
    chn = MIMOChannel(ChnConfig)

    xn = gen.generate_symbols()
    out = chn.process(xn)


    ax1 = plt.subplot(221)
    plt.title("h11")
    fft_plot(chn.coefficients[0,0], ChnConfig.brate)
    plt.grid()
    ax2 = plt.subplot(222, sharey=ax1)
    plt.title("h12")
    fft_plot(chn.coefficients[0,1], ChnConfig.brate)
    plt.grid()
    ax3 = plt.subplot(223, sharex=ax1)
    plt.title("h21")
    fft_plot(chn.coefficients[1,0], ChnConfig.brate)
    plt.grid()
    ax4 = plt.subplot(224, sharey=ax3, sharex=ax1)
    plt.title("h22")
    fft_plot(chn.coefficients[1,1], ChnConfig.brate)
    plt.grid()

    plt.show()
