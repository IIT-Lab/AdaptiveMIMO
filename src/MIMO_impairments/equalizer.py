from configurable import Configurable
# from plotable import Plotable
from slicer import Slicer
import numpy as np


class LMSEqualizer(Configurable):
    def __init__(self, ConfigClass):
        self.taps = 15
        self.dualpol = True
        #Creo matriz de (dualpol+1, taps) i.e cada fila es un vector fila de coeficientes
        self._coefficients = np.zeros((self.dualpol+1, self.taps), dtype=np.complex_)
        self.mu = 1e-4
        self.os_factor = 1
        self.load_config(ConfigClass)

    @property
    def coefficients(self):
        return self._coefficients[:,::-1]

    @coefficients.setter
    def coefficients(self, value):
        self._coefficients = value[:,::-1]

    def filter_once(self, xn):
        return np.inner(self._coefficients, xn)

    def update(self, xn, error):
        self._coefficients += self.mu*np.conj(xn)*error

    def adapt(self, data, reference, **kwargs):


        if isinstance(reference, Slicer):
            diff = lambda x, i: reference.get_symbols(x) - x
        else:
            diff = lambda x, i: reference[:,i].reshape(x.shape) - x

        prefix = int(np.ceil((self.taps-1)/2))
        sufix =  int(np.floor((self.taps-1)/2))



        padded_data = np.zeros((data.shape[0], data.shape[1]+prefix+sufix), dtype=np.complex_)
        padded_data[:,prefix:-sufix] = data

        for i in range(0, data.shape[1], self.os_factor):
            current_data = padded_data[:,i:i+self.taps]
            ye = self.filter_once(current_data)
            coefficients = np.copy(self.coefficients)
            error = diff(ye, i//self.os_factor)
            self.update(current_data, error)
            for key in kwargs:

                kwargs[key].append(locals()[key])


class MIMOEqualizer(LMSEqualizer):
    def __init__(self, ConfigClass):
        super().__init__(ConfigClass)

    @property
    def coefficients(self):
        return self._coefficients[:,:,::-1]

    @coefficients.setter
    def coefficients(self, value):
        self._coefficients = value[:,:,::-1]

    def filter_once(self, xn):
        return np.array([
            [np.tensordot(self._coefficients[0,:], xn, axes=2)],
            [np.tensordot(self._coefficients[1,:], xn, axes=2)],
        ])

    def update(self, xn, error):
        self._coefficients[0,:] += self.mu*np.conj(xn)*error[0]
        self._coefficients[1,:] += self.mu*np.conj(xn)*error[1]

if __name__ == '__main__':
    eq = MIMOEqualizer(None)
    taps = 5
    dualpol = True
    eq.coefficients = np.zeros((dualpol+1,2, taps))
    eq.coefficients[0,0,:] = 0*np.ones(taps)
    eq.coefficients[0,1,:] = 1*np.ones(taps)
    eq.coefficients[1,0,:] = 2*np.ones(taps)
    eq.coefficients[1,1,:] = 3*np.ones(taps)
    xn = np.arange(10).reshape(2,5)

    asd = np.zeros((2,1))
    for i in range(2):
        for j in range(2):
            for k in range(taps):
                asd[i] += eq.coefficients[i,j,k]*xn[j,k]
    print("van los coef")
    print(eq.coefficients)
    print(type(eq.coefficients[0,1][0]))
    print("chau coef")
    print(asd == eq.filter_once(xn))
    print(eq.filter_once(xn))
