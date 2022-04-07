import cygrid
import numpy as np
from read_crafts_spec import read_crafts_spec
from astropy.io import fits
import astropy.units as u

def gridding(speclist,center,imsize=[100,100],cell=1/60,rest_freq=1.420405751768,vmin=-700,vmax=700,p=0,outname='image.fits',overwrite=False):

    c = 299792458.0/1e3 # km/s

    target_header = {
     'NAXIS':  2,
     'NAXIS1': imsize[0],
     'NAXIS2': imsize[1],
     'CTYPE1': 'RA---SIN',
     'CTYPE2': 'DEC--SIN',
     'CUNIT1': 'deg',
     'CUNIT2': 'deg',
     'CDELT1': -cell,
     'CDELT2': cell,
     'CRPIX1': imsize[0] / 2.,
     'CRPIX2': imsize[1] / 2.,
     'CRVAL1': center[0],
     'CRVAL2': center[1],
     }

    print(target_header)
    
    gridder       = cygrid.WcsGrid(target_header)
    FAST_beam     = 3/60  # unit: deg
    kernel_sigma  = FAST_beam/2/np.sqrt(8 * np.log(2))
    kernel_args   = ('gauss1d',(kernel_sigma,),3.*kernel_sigma,kernel_sigma/2.)
    gridder.set_kernel(*kernel_args)
    
    fmin = (1-vmax/c)*rest_freq
    fmax = (1-vmin/c)*rest_freq # unit: GHz

    speclist = np.loadtxt(speclist,dtype=str)
    for specfits in speclist:
        try:
            _,spec_hdu1 = read_crafts_spec(specfits)
        except:
            print('WARNING: BAD FITS FILE!! Please check '+ specfits)
        else:
            print('Gridding '+specfits)
            spec_fmin   = spec_hdu1.data['freq'][0]
            spec_chw    = spec_hdu1.data['chan_bw'][0]
            spec_nch    = spec_hdu1.data['nchan'][0]
            spec_fmax   = spec_fmin+spec_nch*spec_chw
            spec_freq   = np.arange(spec_fmin,spec_fmax,spec_chw)/1e3
            spec_ra     = spec_hdu1.data['OBJ_RA']
            spec_dec    = spec_hdu1.data['OBJ_DEC']
            freq_range  = (spec_freq>fmin) & (spec_freq<fmax) 
            spec_data   = spec_hdu1.data['data'][:,freq_range,p]
            spec_freq   = spec_freq[freq_range]
            gridder.grid(spec_ra, spec_dec, spec_data)
    
    gridded_data = gridder.get_datacube()
    
    hdu = fits.PrimaryHDU(data=gridded_data)
    
    for key in target_header.keys():
        hdu.header[key] = target_header[key]
    
    hdu.header['CTYPE3']  = 'FREQ'
    hdu.header['CUNIT3']  = 'Hz'
    hdu.header['CDELT3']  = spec_chw*1e6
    hdu.header['CRPIX3']  = 0
    hdu.header['CRVAL3']  = spec_freq[0]*1e9
    hdu.header['SPECSYS'] = 'LSRK'
    hdu.header['RESTFRQ'] = rest_freq*1e9 
    hdu.header['BUNIT']   = 'K'

    hdu.writeto(outname,overwrite=overwrite)
