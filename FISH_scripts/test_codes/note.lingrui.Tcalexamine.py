# check the Noise diode spectrum of each beam

# ---------- import packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# ---------- read in Tcal fits
tcalfile = r'CAL.20201014.low.W.fits'
tcalfits = fits.open(tcalfile)

hdu0 = tcalfits[0] # empty hdu
hdu1 = tcalfits[1]

freq = hdu1.data['freq']
tcal = hdu1.data['tcal']

tcal0 = tcal[0,:,:,0]
tcal1 = tcal[0,:,:,1]

beams = tcal.shape[1]

plt.figure()
for i in range(beams):
    plt.step(freq[0],tcal0[i],linestyle='-',where='mid')
    plt.step(freq[0],tcal1[i],linestyle='--',where='mid')
plt.xlabel('Freq')
plt.ylabel('Noise')
plt.ylim(0.9,1.3)
plt.savefig('Noise_diode.pdf',format='pdf')
plt.close()

