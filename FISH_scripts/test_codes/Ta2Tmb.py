import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from read_crafts_spec import read_crafts_spec
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from astropy.utils import iers
#iers.conf.auto_download = False  
iers.conf.auto_max_age = None


def eta_A_curve(zen,beam):
    
    eta_A_file = np.loadtxt('/bb8/lingrui/my_projects/FAST_pipeline/FISH_scripts/eta_A_curve.txt',dtype=str)
    
    i = int(beam-1)
    a = float(eta_A_file[3*i,2])*1e-4
    b = float(eta_A_file[3*i+1,2])*1e-1
    c = float(eta_A_file[3*i+2,2])*1e-2
    d = b + 26.4*(a-c)
    
    if zen <= 26.4:
        return a*zen+b
    if 26.4 < zen <= 40:
        return c*zen+d
    else:
        print('WARNING: BAD ZEN!!! Please Check!!!')

def Ta2Tmb(specfits,beam=1,restfreq=1420.405751768):
    
    spec_hdu0,spec_hdu1 = read_crafts_spec(specfits)
    
    spec_mjd  = Time(np.nanmean(spec_hdu1.data['utobs']),format='mjd')
    spec_ra   = np.nanmean(spec_hdu1.data['obj_ra'])
    spec_dec  = np.nanmean(spec_hdu1.data['obj_dec'])
        
    FAST = EarthLocation(x=-1668557.9202*u.m,y=5506865.964*u.m,z=2744946.0867*u.m)
    
    coord = SkyCoord(ra=spec_ra*u.deg,dec=spec_dec*u.deg,frame='icrs',location=FAST,obstime=spec_mjd)
    
    spec_zen = coord.altaz.zen.value
    
    eta_A  = eta_A_curve(spec_zen,beam=beam)
    wl     = 299792458.0/(restfreq*1e6) # restfreq in MHz
    eta_mb = 0.8899*0.0008436**2*(300/wl)**2*eta_A # The effective diameter of FAST is 300 meter.
                                                   # The beam size at 1420 MHz is ~ 2.9 ' 
                                                   # (~0.0008436 rad)
    
    spec_hdu1.data['data'] = spec_hdu1.data['data']/eta_mb
    
    hduout = fits.HDUList([spec_hdu0,spec_hdu1])
    rootname  = specfits[:-5]
    hduout.writeto(rootname+'_Tmb.fits', overwrite=True)
    hduout.close()

    
    