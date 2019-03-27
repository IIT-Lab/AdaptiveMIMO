import numpy as np
from configurable import Configurable

class SymbolGenerator(Configurable):
    def __init__(self, ConfigClass):
        self.n_tx_symbs = 4e6
        self.M = 2
        self.complex_type = True
        self.dualpol = True
        self.os_factor = 1
        self.load_config(ConfigClass)

    def generate_symbols(self):
        tx_symbols = 2*np.random.randint(
            self.M, size=(self.dualpol+1, int(self.n_tx_symbs))
        )
        tx_symbols -= self.M-1

        if self.complex_type:
            img = 2*np.random.randint(
                self.M, size=(self.dualpol+1, int(self.n_tx_symbs))
            )
            img -= self.M-1
            tx_symbols = tx_symbols + 1j*img

            tx_symbols_os = np.zeros(
                (tx_symbols.shape[0], tx_symbols.shape[1]*self.os_factor),
                dtype=complex
            )
        else:
            tx_symbols_os = np.zeros(
                (tx_symbols.shape[0], tx_symbols.shape[1]*self.os_factor),
            )

        tx_symbols_os[:,0::self.os_factor] = tx_symbols
        return tx_symbols_os
