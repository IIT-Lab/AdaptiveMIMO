import scipy.io as scio

data = scio.loadmat('chn_taps.mat')

h11 = data['h11'].flatten()
h12 = data['h12'].flatten()
h21 = data['h21'].flatten()
h22 = data['h22'].flatten()
