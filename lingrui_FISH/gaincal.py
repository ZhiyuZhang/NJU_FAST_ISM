# Generate the gain table: gain.fits

import numpy as np
import sys
from read_crafts_cal import read_crafts_cal 
import matplotlib.pyplot as plt
from astropy.io import fits

def gaincal(callist,outname='gain',plot=False):
    
    callist = np.loadtxt(callist,dtype=str)
    
    cal_mjd = []
    cal_on  = []
    cal_off = []
    
    for calfile in callist:
    
        hdu0,hdu1,hdu2 = read_crafts_cal(calfile)
    
        cal_mjd.append(hdu2.data['mjd'])
        cal_on.append(hdu2.data['cal'][:,1,:,:])
        cal_off.append(hdu2.data['cal'][:,0,:,:])
    
    cal_mjd = np.concatenate(cal_mjd,axis=0)
    cal_on  = np.concatenate(cal_on,axis=0)
    cal_off = np.concatenate(cal_off,axis=0)
    gain    = cal_on-cal_off
    nsubint = len(cal_mjd)
    
    obsfreq  = hdu0.header['OBSFREQ']                # observing frequency
    obsnchan = hdu0.header['OBSNCHAN']               # observing channel number
    obsbw    = hdu0.header['OBSBW']                  # observing band width
    fmin     = float(obsfreq - obsbw/2.)             # minimum frequency
    fmax     = float(obsfreq + obsbw/2.)             # maximum frequency
    nf       = obsnchan                              # channel number
    df       = hdu1.header['CHAN_BW']                # channel width
    freq     = np.arange(fmin,fmax,df)/1e3
    p        = gain.shape[1]
    
    CAL = fits.Column(name='GAIN', format="%dE" % (nf*p), array=gain, dim=str((nf, p)))
    TIME = fits.Column(name='MJD', format="D", array=cal_mjd)
    cols = fits.ColDefs([TIME, CAL])
    caltab = fits.BinTableHDU.from_columns(cols)
    caltab.writeto(outname + '.fits', overwrite=True)

    # ------------------------- plot the results
    if plot:
        for i in range(gain.shape[1]):
            NoiInj = gain[:,i,:]
            fig = plt.figure()
            gs  = fig.add_gridspec(nrows=2,ncols=2,hspace=0,wspace=0,height_ratios=np.array([2,1]),width_ratios=np.array([2,1]))
            ax00   = fig.add_subplot(gs[0,0])
            ax01   = fig.add_subplot(gs[0,1],sharey=ax00)
            ax10   = fig.add_subplot(gs[1,0],sharex=ax00)
            Noimap = ax00.imshow(NoiInj,cmap='jet',aspect='auto',vmin=0,vmax=1)
            ax10.step(np.arange(0,nf),np.average(NoiInj,axis=0),where='mid')
            ax01.step(np.average(NoiInj,axis=1),np.arange(0,nsubint),where='mid')
            ax00.tick_params(labelbottom=False)
            ax01.tick_params(labelleft=False)
            ax00.set_ylabel('Subintegration')
            ax01.set_xlabel('ON-OFF')
            ax10.set_xlabel('Channel')#,fontweight='bold')
            ax10.set_ylabel('ON-OFF')
            
            plt.savefig('gaincurveVStime_{}.pdf'.format(i))
            plt.close()
    return "## End of Task GAINCAL ##"

if __name__ == '__main__':
    callist = sys.argv[1]
    plot = sys.argv[2]
    gaincal(callist,plot)
