from read_crafts_spec import read_crafts_spec
from astropy.io import fits
#def fast2class(specfits):
specfits = 'Calibrated/Dec+2150_06_05_arcdrift-M01_N_0001_210204_135403to210204_140055_calibrated_icrs.fits'
spec_hdu0,spec_hdu1 = read_crafts_spec(specfits)

spec_hdu1_new = spec_hdu1

spec_hdu1_new.header['MAXIS']  = 4
spec_hdu1_new.header['MAXIS1'] = spec_hdu1.data['nchan'][0]
spec_hdu1_new.header['MAXIS1'] = spec_hdu1.data['nchan'][0]
