import numpy as np
#import pyfits
import astropy.io.fits as pyfits
from astropy.time import Time as ATime
from datetime import timedelta
import os, sys
from decimal import Decimal
import time

def extract_crafts_cal(filename):
    
    Nbin = 2 
    Tconst = 0.201326592
    
    hdulist = pyfits.open(filename)
    hdu0 = hdulist[0]
    hdu1 = hdulist[1]
    
    #header = hdu0.header
    obsfreq = hdu0.header['OBSFREQ']
    obsnchan = hdu0.header['OBSNCHAN']
    obsbw = hdu0.header['OBSBW']
    fmin = float(obsfreq - obsbw/2.)
    fmax = float(obsfreq + obsbw/2.)
    nf = obsnchan
    df = hdu1.header['CHAN_BW']
    tsamp = hdu1.header['TBIN']
    nsubint = hdu1.header['NAXIS2']
    samppersubint = int(hdu1.header['NSBLK'])
    nsamp = nsubint * samppersubint
    sourename = hdu0.header['SRC_NAME']
    #ra = hdu0.header['RA']
    #dec = hdu0.header['DEC']
    nbits = hdu0.header['BITPIX']
    header = dict(hdu0.header + hdu1.header)
    duration = tsamp * nsamp
    #Nsub = nsubint
    #Tsub = tsamp * samppersubint
    
    samppersubint = int(Tconst / tsamp)
    Nsub = int(nsamp / samppersubint)
    Tsubint = Tconst
    
    obsmjd = int(hdu0.header['STT_IMJD'])
    #tstart = Decimal(hdu0.header['STT_SMJD']) + Decimal(hdu0.header['STT_OFFS'])
    tstart = float(hdu0.header['STT_SMJD']) + float(hdu0.header['STT_OFFS'])
    
    #mjd = obsmjd + tstart/Decimal(24*3600)
    mjd = obsmjd + tstart/86400 ## start mjd
    
    header['Tsubint'] = Tsubint
    header['Duration'] = duration
    header['obsmjd'] = str(mjd)
    obsdate = ATime(obsmjd, format='mjd').datetime
    #obsdate.format = 'datetime'
    #print obstime.datetime + timedelta( seconds = float(tstart) )
    #print  timedelta( seconds = float(tstart) )
    obstime = timedelta( seconds = float(tstart) )
    obsutc = obsdate + obstime
    header['obsutc'] = obsutc
    
    print('tsamp', tsamp)
    print('mjd', mjd)
    print('obsutc', header['obsutc'].strftime("%y-%m-%d %H:%M:%S") )
    print('Duration:', duration)
    print('Tsubint:', Tsubint)
    #sys.exit(0)
    
    data = np.copy(hdu1.data['data'][:,:,:,:,:])
    
    print('data shape before fold', data.shape)
    l,m,p,n,o = data.shape
    Npol = p
    
    if Npol > 2:
        dataauto = data[:,:,:2,:,:]
        datacros = data[:,:,2:,:,:].astype(np.int8)
        data0 = dataauto.reshape(Nsub,-1,Nbin,2,n).mean(axis=1, dtype=np.float32) # group the 128 subintegrations two-by-two; group the 2*1024 sampling two-by-two; average in each subintegration group
        data1 = datacros.reshape(Nsub,-1,Nbin,2,n).mean(axis=1, dtype=np.float32)
        data = np.concatenate((data0, data1), axis=2)
    else:
        data = data.reshape(Nsub,-1,Nbin,p,n)
        data = data.mean(axis=1, dtype=np.float32)
    #CALON = data[:,1,:,:]
    #CALOFF = data[:,0,:,:]
    #CALDATA = np.stack((CALON, CALOFF))
    Tarray = mjd + (np.arange(Nsub) + 0.5) * Tconst / 86400 # center time of each sampling
    
    #now = time.time()
    
    delnames = ['DATA', 'DAT_OFFS', 'DAT_SCL', 'DAT_WTS', 'DAT_FREQ']
    datadict = {}
    for name in hdu1.data.names:
        if name not in delnames:
            datadict[name] = hdu1.data[name]
    
    hdu1.columns.del_col('DATA')
    hdu1.columns.del_col('DAT_OFFS')
    hdu1.columns.del_col('DAT_SCL')
    hdu1.columns.del_col('DAT_WTS')
    hdu1.columns.del_col('DAT_FREQ')
    pyfits.BinTableHDU.update(hdu1)
    for name in datadict:
        hdu1.data[name] = datadict[name]
    
    #print 'took %f second to form mid block' % ( time.time() - now)
    
    CAL = pyfits.Column(name='CAL', format="%dE" % (Nbin*p*n), array=data, dim=str((n, p, Nbin)))
    TIME = pyfits.Column(name='MJD', format="D", array=Tarray)
    cols = pyfits.ColDefs([TIME, CAL])
    caltab = pyfits.BinTableHDU.from_columns(cols)
    hdulist.append(caltab)
    
    
    
    rootname = filename[:-5]
    decstr = rootname.split('_')[0]
    beamnum = rootname.split('-M')[-1]
    newname = decstr + '-M' + beamnum + "_" + obsutc.strftime("%y%m%d") + "_" + obsutc.strftime("%H%M%S")
    
    #np.savez(Savepath + rootname + '.npz', header=header, Nbin=Nbin, Nsub=Nsub, Npol=Npol, data = data.astype(np.float33))
    #np.savez(newname + '.npz', header=header, Nbin=Nbin, Nsub=Nsub, Npol=Npol, dimensions=('Sub integration', 'On and Off', 'Polarization', 'Spectral Channels'), data = data.astype(np.float))
    hdulist.writeto(newname + '.fits.gz', overwrite=True)


if __name__ == '__main__':
    filename = sys.argv[1]
    extract_crafts_cal(filename)
