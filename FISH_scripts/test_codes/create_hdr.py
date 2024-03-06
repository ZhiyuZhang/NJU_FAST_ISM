from argparse import _MutuallyExclusiveGroup
from read_crafts_spec import read_crafts_spec

specfits = r'Calibrated/Dec+2150_06_05_arcdrift-M01_N_0001_210204_135403to210204_140055_calibrated.fits'
_,spec_hdu1 = read_crafts_spec(specfits)
imsize = [60, 80]
cell   = 1/60
center = [94.2850057, 22.5402125] 
spec_fmin = spec_hdu1.data['freq'][0]
spec_chw  = spec_hdu1.data['chan_bw'][0]
spec_nch  = spec_hdu1.data['nchan'][0]

target_header = {
     'NAXIS': 3,
     'NAXIS1': imsize[0],
     'NAXIS2': imsize[1],
     'NAXIS3': spec_nch,
     'CTYPE1': 'RA---SIN',
     'CTYPE2': 'DEC--SIN',
     'CTYPE3': 'FREQ',
     'CUNIT1': 'deg',
     'CUNIT2': 'deg',
     'CUNIT3': 'Hz',
     'CDELT1': -cell,
     'CDELT2': cell,
     'CDELT3': spec_chw*1e6,
     'CRPIX1': imsize[0] / 2.,
     'CRPIX2': imsize[1] / 2.,
     'CRPIX3': 0,
     'CRVAL1': center[0],
     'CRVAL2': center[1],
     'CRVAL3': spec_fmin*1e6,
     }


speclist = 'mygridlist_all'
gridding(speclist,center=[94.2850057, 22.5402125],imsize=[90,90],cell=1/60,rest_freq=1.420405751768,vmin=-500,vmax=500,p=0,outname='image.fits',overwrite=True)

import numpy as np
ra = []
dec = []
speclist = np.loadtxt(speclist,dtype=str)
for specfits in speclist:
     _,spec_hdu1 = read_crafts_spec(specfits)
     ra.append(spec_hdu1.data['OBJ_RA'])
     dec.append(spec_hdu1.data['OBJ_DEC'])

for i in range(len(ra)):
     plt.plot(ra[i],dec[i],'o')

