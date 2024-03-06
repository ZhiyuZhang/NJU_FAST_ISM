# transform the (x,y,z) coordinates of the feed to celestial frame

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from read_crafts_spec import read_crafts_spec
from read_crafts_cal import read_crafts_cal
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from astropy.utils import iers
#iers.conf.auto_download = False  
iers.conf.auto_max_age = None  

def interp(x,y,xnew,kind='linear',axis=0):
    interp_func = interp1d(x=x,y=y,kind=kind,axis=axis,fill_value='extrapolate')
    ynew = interp_func(xnew)
    return ynew

def xyz2celestial(specfits,xyztable,beam=1,plot=False,psrflag=False):

    if not psrflag:
        spec_hdu0,spec_hdu1 = read_crafts_spec(specfits)
        spec_mjd  = spec_hdu1.data['utobs']
    else:
        spec_hdu0,spec_hdu1,spec_hdu2 = read_crafts_cal(specfits)
        spec_mjd = spec_hdu2.data['MJD']
    
    #FAST = EarthLocation(lon=106.85645657571428*u.deg,lat=25.6534387*u.deg,height=1138.72178*u.m)
    FAST = EarthLocation(x=-1668557.9202*u.m,y=5506865.964*u.m,z=2744946.0867*u.m)
    
    cabin_time = np.array(xyztable['SysTime'].values,dtype=str)
    cabin_utc  = Time(cabin_time,format='iso')-8*u.hr
    cabin_mjd  = cabin_utc.mjd
    range      = (cabin_mjd>=spec_mjd[0]) & (cabin_mjd<=spec_mjd[-1])
    cabin_mjd  = cabin_mjd[range]
    
    ## East-North-Sky System
    SDP_E = xyztable['SDP_PhaPos_X'].values[range]
    SDP_N = xyztable['SDP_PhaPos_Y'].values[range]
    SDP_S = xyztable['SDP_PhaPos_Z'].values[range]
    
    SDP_E_interp = interp(cabin_mjd,SDP_E,spec_mjd)
    SDP_N_interp = interp(cabin_mjd,SDP_N,spec_mjd)
    SDP_S_interp = interp(cabin_mjd,SDP_S,spec_mjd)
    
    spec_time  = Time(spec_mjd,format='mjd')
    coord      = SkyCoord(x=-SDP_N_interp,y=-SDP_E_interp,z=-SDP_S_interp,representation_type='cartesian',frame='altaz',location=FAST,obstime=spec_time)
    coord.representation_type = 'spherical'

    MBeams    = {'01':(0     , 0)     , '02':(5.74    , 0.00811) , '03':(2.88 , -4.97)  , '04':(-2.86 , -4.98) , '05':(-5.74 , -0.0127) ,
                 '06':(-2.88 , 4.97)  , '07':(2.86    , 4.98)    , '08':(11.5 , 0.0116) , '09':(8.61  , -4.96) , '10':(5.75  , -9.93)   ,
                 '11':(0.018 , -9.94) , '12':(-5.71   , -9.96)   , '13':(-8.6 , -4.99)  , '14':(-11.5 , -0.03) , '15':(-8.63 , 4.95)    ,
                 '16':(-5.77 , 9.93)  , '17':(-0.0181 , 9.94)    , '18':(5.73 , 9.95)   , '19':(8.61  , 4.98)}
    SDP_angle = np.unique(xyztable['SDP_AngleM'].values[range])[0]

    beam    = '{}'.format(str(beam).rjust(2,'0'))
    pos     = MBeams[beam]
    rot     = np.array([[np.cos(SDP_angle),-np.sin(SDP_angle)],[np.sin(SDP_angle),np.cos(SDP_angle)]])
    pos_rot = np.dot(rot,pos)

    obj_ra  = coord.icrs.ra.value + pos_rot[0]/np.cos(coord.icrs.dec.value*u.deg)/60
    obj_dec = coord.icrs.dec.value + pos_rot[1]/60
    
    n = len(obj_ra)
    
    if not psrflag:
        spec_hdu1.data['OBJ_RA']  = obj_ra
        spec_hdu1.data['OBJ_DEC'] = obj_dec
        hduout = fits.HDUList([spec_hdu0,spec_hdu1])
    else:
        col_ra  = fits.Column(name='OBJ_RA',format="1D", array=obj_ra)
        col_dec = fits.Column(name='OBJ_DEC',format="1D", array=obj_dec)
        spec_hdu2.columns.add_col(col_ra)
        spec_hdu2.columns.add_col(col_dec)
        hduout = fits.HDUList([spec_hdu0,spec_hdu1,spec_hdu2])
    
    rootname  = specfits[:-5]
    hduout.writeto(rootname+'_icrs.fits', overwrite=True)
    hduout.close()

    if plot:
        fig = plt.figure()
        plt.plot(spec_hdu1.data['OBJ_RA'],spec_hdu1.data['OBJ_DEC'],color='black')
        plt.xlabel('R.A. [ICRS]')
        plt.ylabel('Dec. [ICRS]')
        plt.savefig(rootname+'_scan.pdf',format='pdf',bbox_inches='tight')
        plt.close()

#def xyztoradec(xyz):
#    R=np.sqrt(xyz[0]**2+xyz[1]**2+xyz[2]**2)
#    z=-xyz[2]
#    y=-xyz[0]
#    x=-xyz[1]
#    az0 = 0.0
#    if (np.abs(x) < 1.0e-8):                 # If X is a very small number
#       if (y>0):
#          az0 = 90.0                         # abs(X) <1E-8; Y >0;  AZ = 90
#       else:
#          az0 = 270.0                        # abs(X) <1E-8; Y <=0; AZ = 270
#    else:
#       tempaz0 = np.arctan(y/x)*180/np.pi      # X >=1E-8;  arctan(y/x) to degree
#       if (x>0 and y>0):
#          az0 = tempaz0                      # X > 0; Y > 0; AZ = tempaz0
#       if ((x>0 and y<0) or (x>0 and y==0)):
#           az0 = tempaz0+360.0               # X > 0; Y <= 0;   AZ = tempaz0+360
#       if (x<0 and y>0):
#           az0 = tempaz0+180.0               # X < 0; Y > 0;    AZ = tempaz0+180
#       if ((x<0 and y<0) or (x<0 and y==0)):
#           az0 = tempaz0+180.0               # X < 0; Y <= 0;   AZ = tempaz0+180
#    el0 = np.arcsin(z/R)*180/np.pi
#    return az0,el0

