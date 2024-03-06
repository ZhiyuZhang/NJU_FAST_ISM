# apply the TSYS_PSR

import numpy as np
import sys
from read_crafts_cal import read_crafts_cal
from read_crafts_spec import read_crafts_spec
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time as ATime
from scipy.interpolate import interp1d
from scipy import ndimage

def interp(x,y,xnew,kind='linear',axis=0,smooth=0,smoothmode='nearest'):
    if smooth != 0:
        y = ndimage.uniform_filter1d(y,smooth,axis=axis,mode=smoothmode)
    interp_func = interp1d(x=x,y=y,kind=kind,axis=axis,fill_value='extrapolate')
    ynew = interp_func(xnew)
    return ynew

def applycal(specfits,tsysfits,restfreq=1.420405751768,plot=False):

    print('Matching the Tsys of {} and {}'.format(specfits,tsysfits))
    
    spec_hdu0,spec_hdu1 = read_crafts_spec(specfits)
    spec_mjd  = spec_hdu1.data['UTOBS']
    spec_data = spec_hdu1.data['DATA']
    spec_fmin = spec_hdu1.data['freq'][0]
    spec_chw  = spec_hdu1.data['chan_bw'][0]
    spec_nch  = spec_hdu1.data['nchan'][0]
    spec_fmax = spec_fmin+spec_nch*spec_chw
    spec_freq = np.arange(spec_fmin,spec_fmax,spec_chw)/1e3
    
    cal_hdu0,cal_hdu1,cal_hdu2 = read_crafts_cal(tsysfits)
    cal_mjd  = cal_hdu2.data['MJD']
    tpsr     = cal_hdu2.data['TSYS_PSR']
    obsfreq  = cal_hdu0.header['OBSFREQ']  # observing frequency
    obsnchan = cal_hdu0.header['OBSNCHAN'] # observing channel number
    obsbw    = cal_hdu0.header['OBSBW']    # observing band width
    fmin     = float(obsfreq - obsbw/2.)   # minimum frequency
    fmax     = float(obsfreq + obsbw/2.)   # maximum frequency
    nf       = obsnchan                    # channel number
    df       = cal_hdu1.header['CHAN_BW']  # channel width
    cal_freq = np.arange(fmin,fmax,df)/1e3 # frequency array (GHz)
    
    tpsr_interp = interp(cal_mjd,tpsr,spec_mjd,kind='linear',axis=0)
    freq_range  = np.array([-700,700])/3e5*restfreq+restfreq  # -700 km/s to 700 km/s; 
                                                              # ignore this freq_range in matching
    spec_median = np.nanmedian(spec_data[:,(spec_freq<freq_range[0]) | (spec_freq>freq_range[1]),:],axis=1)
    tpsr_median = np.nanmedian(tpsr_interp[:,:,((cal_freq>spec_freq[0]) & (cal_freq<freq_range[0])) | ((cal_freq>freq_range[1]) & (cal_freq<spec_freq[-1]))],axis=2)

    factor      = spec_median/tpsr_median
    factor      = np.transpose(np.tile(factor,(1,1,1)),axes = (1,0,2))
    spec_data_T = spec_data/factor
    
    del spec_data
    #T_A    = fits.Column(name='T_A', format="%dE" % (nf*p), array=spec_data_T,dim=str((p,nf)))
    #TIME   = fits.Column(name='MJD', format="D", array=spec_mjd)
    #Tatab  = fits.BinTableHDU.from_columns(fits.ColDefs([T_A,TIME]))
    #keys = ['TELESCOP','BANDWID','INSTRUME','CTYPE1','CTYPE2']
    #for k in keys:
    #    Tatab.header[k] = hdu1.header[k]
    #Tatab.header['FREQ'] = spec_fmin*1e6
    #Tatab.header['CHAN_BW'] = spec_chw*1e6
    #Tatab.header['NCHAN'] = spec_nch
    spec_hdu1.data['DATA'] = spec_data_T

    hduout = fits.HDUList([spec_hdu0,spec_hdu1])
    rootname  = specfits[:-5].split('/')[-1]
    spec_date = ATime(spec_mjd,format='mjd').datetime
    newname   = 'Calibrated/' + rootname + '_' + spec_date[0].strftime("%y%m%d") + "_" + spec_date[0].strftime("%H%M%S") + 'to' \
              + spec_date[-1].strftime("%y%m%d") + "_" + spec_date[-1].strftime("%H%M%S")
    hduout.writeto(newname+'_calibrated.fits', overwrite=True)
    hduout.close()

    if plot:
            for i in range(2):
                SPEC = spec_data_T[:,(spec_freq>freq_range[0]) & (spec_freq<freq_range[1]),i]
                nsubint,nf = SPEC.shape
                fig  = plt.figure()
                gs   = fig.add_gridspec(nrows=2,ncols=2,hspace=0,wspace=0,height_ratios=np.array([2,1]),width_ratios=np.array([2,1]))
                ax00   = fig.add_subplot(gs[0,0])
                ax01   = fig.add_subplot(gs[0,1],sharey=ax00)
                ax10   = fig.add_subplot(gs[1,0],sharex=ax00)
                Noimap = ax00.imshow(SPEC,cmap='jet',aspect='auto')
                ax10.step(np.arange(0,nf),np.nanmedian(SPEC,axis=0),where='mid')
                ax01.step(np.nanmedian(SPEC,axis=1),np.arange(0,nsubint),where='mid')
                ax00.tick_params(labelbottom=False)
                ax01.tick_params(labelleft=False)
                ax00.set_ylabel('Subintegration')
                ax01.set_xlabel(r'$T_{\rm A}$')
                ax10.set_xlabel('Channel')#,fontweight='bold')
                ax10.set_ylabel(r'$T_{\rm A}$')
                cax = fig.add_axes([ax00.get_position().x0,ax00.get_position().y1,ax00.get_position().width,0.02])
                cb  = fig.colorbar(mappable=Noimap,cax=cax,pad=0,orientation='horizontal')
                cb.ax.tick_params(labelsize=15,direction='in',labeltop=True,labelbottom=False)
                cax.set_title(r'$T_{\rm A}$',fontsize=15)
                plt.savefig('Calibrated/'+rootname+'_p{}.pdf'.format(i),format='pdf',bbox_inches='tight')
                plt.close()

