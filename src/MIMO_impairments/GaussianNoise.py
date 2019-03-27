# GaussianNoise.py
import numpy as np
import matplotlib.pyplot as plt

def error_calculation(bit_length, noise_amp):
  #uniformly distributed sequence to generate singal
  b = np.random.uniform(-1, 1, bit_length)

  #binary signal generated from 'b'
  signal = np.zeros((bit_length),float)
  for i in range(len(b)):
    if b[i] < 0:
      signal[i]=-1
    else:
      signal[i]=1

  #Gaussian Noise
  noise = np.random.randn(bit_length)

  #recevied signal
  rx_signal = signal + noise_amp*noise

  detected_signal = np.zeros((bit_length),float)
  for i in range(len(b)):
    if rx_signal[i] < 0:
      detected_signal[i]=-1
    else:
      detected_signal[i]=1

  error_matrix = abs((detected_signal - signal)/2)
  error=error_matrix.sum()
  return error