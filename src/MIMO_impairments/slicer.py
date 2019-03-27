from configurable import Configurable
import numpy as np

class Slicer(Configurable):
    def __init__(self, ConfigClass):
        self.os_factor = 1
        self.M = 2

    def get_symbols(self, yn):
        y_ss = yn[::self.os_factor]
        detected_symbols = 2*np.floor(y_ss/2)+1
        return np.clip(detected_symbols, -self.M*2-1, self.M*2+1)
