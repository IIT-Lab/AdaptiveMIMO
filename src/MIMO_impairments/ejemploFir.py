def fixing_append(vector):
    fixed_vector = []
    for ptr in range(len(vector)):
        fixed_vector.append(vector[ptr].fValue)
    return fixed_vector

######################################################################

import math
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
from tool._fixedInt import *
from tool.r_cosine import *
from mpl_toolkits.mplot3d import Axes3D

######################################################################
# Parametros generales
######################################################################
Tbaud     = 1.0/1024000.0
os_tx     = 4.0
beta_tx   = 0.5
Nbauds_tx = 8.0

######################################################################
# Filtro SRRC 
######################################################################

# Coeficientes Filtro transmisor
(t,rrc_tx) = rrcosine(beta_tx, Tbaud, os_tx, Nbauds_tx)                            # Filtro Tx en punto flotante
rrc_tx_fx0  = fixing_append(arrayFixedInt(12,10,rrc_tx,'S','trunc','saturate'))    # Filtro Tx en punto fijo (12,10)
rrc_tx_fx1  = fixing_append(arrayFixedInt(8 ,6 ,rrc_tx,'S','trunc','saturate'))    # Filtro Tx en punto fijo (8,6)
rrc_tx_fx2  = fixing_append(arrayFixedInt(8 ,7 ,rrc_tx,'S','round','saturate'))    # Filtro Tx en punto fijo (8,7)
rrc_tx_fx3  = fixing_append(arrayFixedInt(8 ,8 ,rrc_tx,'S','round','saturate'))    # Filtro Tx en punto fijo (8,8)


fig = plt.figure()
plt.title('Digital filter frequency response')
plt.subplot(211)
w,h = scipy.signal.freqz(rrc_tx)
plt.plot(w, 20 * np.log10(abs(h)),label='Float')
w,h = scipy.signal.freqz(rrc_tx_fx0)
plt.plot(w, 20 * np.log10(abs(h)),label='S(12.10)')
w,h = scipy.signal.freqz(rrc_tx_fx1)
plt.plot(w, 20 * np.log10(abs(h)),label='S(8.6)')
w,h = scipy.signal.freqz(rrc_tx_fx2)
plt.plot(w, 20 * np.log10(abs(h)),label='S(8.7)')
w,h = scipy.signal.freqz(rrc_tx_fx3)
plt.plot(w, 20 * np.log10(abs(h)),label='S(8.8)')
plt.legend(prop={'size': 8})
plt.ylabel('Amplitude [dB]', color='b')
plt.xlabel('Frequency [rad/sample]')
plt.grid()
plt.axis('tight')
plt.ylim(-60,20)


plt.subplot(212)
plt.plot(rrc_tx,'-x',label='Float')
plt.plot(rrc_tx_fx0,'-x',label='S(12.10)')
plt.plot(rrc_tx_fx1,'-x',label='S(8.6)')
plt.plot(rrc_tx_fx2,'-x',label='S(8.7)')
plt.plot(rrc_tx_fx3,'-x',label='S(8.8)')
plt.legend(prop={'size': 8})
plt.ylabel('Amplitude', color='b')
plt.xlabel('Samples')
plt.grid()
plt.axis('tight')
plt.ylim(-0.2,1.2)

plt.savefig('freqandimpulse.eps')

plt.show(block=False)
raw_input('Press Enter to Continue')
plt.close()

rrc_tx_fx0 = arrayFixedInt(12,10,rrc_tx,'S','trunc','saturate')
for ptr in range(len(rrc_tx_fx0)):
    print rrc_tx_fx0[ptr].intvalue

    
# yk_i_fix[ptr].intvalue
# sin_fix._setValue(np.sin(-buff_vec_r[5]))
# DeFixedInt(6,5,'S','trunc','saturate')
