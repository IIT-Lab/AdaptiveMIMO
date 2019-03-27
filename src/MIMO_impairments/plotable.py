import numpy as np
from scipy.fftpack import fft, fftfreq, fftshift
import matplotlib.pyplot as plt

class Plotable:
    def plot(PlotConfig, **kwargs_decorator):
        def decorator_wrapper(func):
            def wrapper(*args, **kwargs):
                result = func(*args, **kwargs, **kwargs_decorator)
                for key, config in zip(kwargs_decorator, PlotConfig):
                    for key in config:
                        getattr(plt, key)(*config[key])
                return result
            return wrapper
        return decorator_wrapper

def fft_plot(y, Fs, db=False, ax=False):
    # number of signal points
    # sample spacing
    ax = ax or plt.gca()
    N = len(y)
    T = 1.0 / Fs
    yf = fft(y)
    xf = fftfreq(N, T)
    xf = fftshift(xf)
    yplot = 1.0/N * np.abs(fftshift(yf))
    if db:
        yplot = 20*np.log10(yplot)
    return ax.plot(xf, yplot)
