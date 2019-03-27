import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy.special import erfc
from tool import _fixedInt as fix
import funciones

NBT = 14
NBF = 12



class LMS_equalizer():
    def __init__(self, n_taps, step):
        #Param: cantidad de taps del filtro
        self.n_taps = int(n_taps)
        self.step   = step

        #Atributo: vector de coeficientes, inicialmente un impulso.
        self.coefficients = np.zeros(self.n_taps, dtype = 'float64')
        self.coefficients[int((n_taps-1)/2)] = 1
        self.c_history = list()
        self.mse = list()

    def filter_once(self, xn):
        return np.array(np.dot(xn, self.coefficients), dtype = 'float64')

    def update(self, xn, error):
        self.coefficients += self.step*error*xn
        self.c_history.append(self.coefficients)

    def adapt(self, input_signal):

        output_eq = []
        error_out = []
        # detected  = []
        # for step in steps:


        for k in range(0, len(input_signal) - self.n_taps, 1):

            # Obtengo muestras
            sample_buffer = np.array(input_signal[k:self.n_taps + k], dtype='float64')


            # Swap para la convolucion
            sample_buffer = sample_buffer[::-1]

            #Output equalizer en k
            yk = self.filter_once(sample_buffer)

            # Computo error
            ek = np.sign(yk) - yk
            self.mse.append(float(ek ** 2))

            #Actualizo valor de coeficientes
            self.update(sample_buffer, ek)

            output_eq.append(yk)
            error_out.append(ek)


        error_out = np.array(error_out)


        return np.array(output_eq), np.array(error_out)






class LMS_equalizer_fixed():

    def __init__(self, n_taps, step):
        #Param: cantidad de taps del filtro
        self.n_taps = int(n_taps)
        self.step   = step

        #Atributo: vector de coeficientes, inicialmente un impulso.
        self.coefficients = np.zeros(self.n_taps, dtype = 'float64')
        self.coefficients[int((n_taps-1)/2)] = 1
        self.coefficients = np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,self.coefficients,'S','round','saturate'),0))
        self.c_history = list()
        self.mse = list()

    def filter_once(self, xn):
        return np.dot(xn, self.coefficients)
        #return np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,np.dot(xn, self.coefficients),'S','round','saturate'),0))


    def update(self, xn, error):
        self.coefficients += self.step*error*xn
        self.coefficients = np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,self.coefficients,'S','round','saturate'),0))
        self.c_history.append(self.coefficients)

    def adapt(self, input_signal):

        output_eq = []
        error_out = []
        # detected  = []
        # for step in steps:


        for k in range(0, len(input_signal) - self.n_taps, 1):

            # Obtengo muestras
            sample_buffer = np.array(input_signal[k:self.n_taps + k])

            # Swap para la convolucion
            sample_buffer = sample_buffer[::-1]

            #Output equalizer en k
            yk = self.filter_once(sample_buffer)

            # Computo error
            ek = np.sign(yk) - yk
            self.mse.append(float(ek ** 2))

            #Actualizo valor de coeficientes
            self.update(sample_buffer, ek)

            output_eq.append(yk)
            error_out.append(ek)
            # if k ==25:
            #     print('Ek FLOAT ES', ek)


        output_eq = np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,output_eq,'S','round','saturate'),0))
        error_out = np.array(funciones.fix_append(fix.arrayFixedInt(NBT,NBF,error_out,'S','round','saturate'),0))
        #print('El error fixed es', error_out[25])


        return output_eq, error_out