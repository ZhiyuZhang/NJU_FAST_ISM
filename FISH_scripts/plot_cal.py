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

# ------------------------- plot on-off and resampled noise spectra

def plot_Tcal(Tcalfits='Tcal/20201014/CAL.20201014.low.W.fits',beam=1):
    
    _,noisehdu1 = read_crafts_spec(Tcalfits)
    noise_data  = np.transpose(noisehdu1.data['tcal'][0,beam-1,:,:],axes=(1,0))
    print('There are {} polarizations in the Tcal noise.'.format(noise_data.shape[0]))
    # if noise_data.shape[0]<p: 
    # # add NaN array when the length of the polarization axis of the noise data is smaller than 4.
    #     noise_data = np.concatenate((noise_data,np.nan*np.ones(shape=(p-noise_data.shape[0],noise_data.shape[1]))),axis=0)
    noise_freq  = noisehdu1.data['freq'][0]/1e3 # MHz to GHz

            
    fig  = plt.figure()
    gs   = fig.add_gridspec(nrows=2,ncols=1,hspace=0,height_ratios=np.array([1,1]))
    ax00 = fig.add_subplot(gs[0,0])
    ax10 = fig.add_subplot(gs[1,0],sharex=ax00)
    ax00.step(noise_freq,noise_data[0],where='mid',label='AA')
    #ax00.step(cal_freq,noise_data_interp[0],linestyle='--',where='mid',label='AA resampled')
    ax10.step(noise_freq,noise_data[1],where='mid',label='BB')
    #ax10.step(cal_freq,noise_data_interp[1],linestyle='--',where='mid',label='BB resampled')
    ax00.set_ylim(1.0,1.3)
    ax10.set_ylim(1.0,1.3)
    #ax00.set_xlim(1.40,1.44)
    ax00.legend()
    ax10.legend()
    ax10.set_xlabel('FREQ [GHz]')
    ax10.set_ylabel('T_CAL [K]')
    ax00.set_title(Tcalfits,fontsize=15,y=0.98)
    plt.savefig('Tcal.pdf',format='pdf',bbox_inches='tight')
    plt.close()
    
    return

def plot_Tpsr(Tpsrfits):
    
    tsys_hdu0,tsys_hdu1,tsys_hdu2 = read_crafts_cal(Tpsrfits)
    
    Tpsr      = tsys_hdu2.data['TSYS_PSR']
    tsys_mjd  = tsys_hdu2.data['MJD']
    obsfreq   = tsys_hdu0.header['OBSFREQ']  # observing frequency
    obsnchan  = tsys_hdu0.header['OBSNCHAN'] # observing channel number
    obsbw     = tsys_hdu0.header['OBSBW']    # observing band width
    fmin      = float(obsfreq - obsbw/2.)   # minimum frequency
    fmax      = float(obsfreq + obsbw/2.)   # maximum frequency
    nf        = obsnchan                    # channel number
    df        = tsys_hdu1.header['CHAN_BW']  # channel width
    cal_freq  = np.arange(fmin,fmax,df)/1e3 # frequency array (GHz)
    nsubint   = len(tsys_mjd)

    for i in range(2):
        T    = Tpsr[:,i,:]
        fig  = plt.figure()
        gs   = fig.add_gridspec(nrows=2,ncols=2,hspace=0,wspace=0,height_ratios=np.array([2,1]),width_ratios=np.array([2,1]))
        ax00 = fig.add_subplot(gs[0,0])
        ax01 = fig.add_subplot(gs[0,1],sharey=ax00)
        ax10 = fig.add_subplot(gs[1,0],sharex=ax00)
        Tmap = ax00.imshow(T,cmap='jet',aspect='auto',vmin=25,vmax=50)
        ax10.step(np.arange(0,nf),np.nanmedian(T,axis=0),where='mid')        
        ax01.step(np.nanmedian(T,axis=1),np.arange(0,nsubint),where='mid')
        ax00.tick_params(labelbottom=False)
        ax01.tick_params(labelleft=False)
        ax00.set_ylabel('Subintegration')
        ax01.set_xlabel('Tpsr')
        ax10.set_xlabel('Channel')#,fontweight='bold')
        ax10.set_ylabel('Tpsr')
        ax10.set_ylim(10,40)
        ax00.set_title(Tpsrfits,fontsize=5)
        cax = fig.add_axes([ax00.get_position().x0,ax00.get_position().y1,ax00.get_position().width,0.02])
        cb  = fig.colorbar(mappable=Tmap,cax=cax,pad=0,orientation='horizontal')
        cb.ax.tick_params(labelsize=15,direction='in',labeltop=True,labelbottom=False)
        cax.set_title(r'$T_{\rm A}$ [K]',fontsize=15)
        plt.savefig('Tpsr2D_p{}.pdf'.format(i),format='pdf',bbox_inches='tight')
        plt.close()

    return  


def plot_off(calfits):
    
    cal_hdu0,cal_hdu1,cal_hdu2 = read_crafts_cal(calfits)
    
    cal_off   = cal_hdu2.data['CAL'][:,0,:,:]
    cal_mjd   = cal_hdu2.data['MJD']
    obsfreq   = cal_hdu0.header['OBSFREQ']  # observing frequency
    obsnchan  = cal_hdu0.header['OBSNCHAN'] # observing channel number
    obsbw     = cal_hdu0.header['OBSBW']    # observing band width
    fmin      = float(obsfreq - obsbw/2.)   # minimum frequency
    fmax      = float(obsfreq + obsbw/2.)   # maximum frequency
    nf        = obsnchan                    # channel number
    df        = cal_hdu1.header['CHAN_BW']  # channel width
    cal_freq  = np.arange(fmin,fmax,df)/1e3 # frequency array (GHz)
    nsubint   = len(cal_mjd)

    for i in range(2):
        T    = cal_off[:,i,:]
        fig  = plt.figure()
        gs   = fig.add_gridspec(nrows=2,ncols=2,hspace=0,wspace=0,height_ratios=np.array([2,1]),width_ratios=np.array([2,1]))
        ax00 = fig.add_subplot(gs[0,0])
        ax01 = fig.add_subplot(gs[0,1],sharey=ax00)
        ax10 = fig.add_subplot(gs[1,0],sharex=ax00)
        Tmap = ax00.imshow(T,cmap='jet',aspect='auto',vmin=25,vmax=50)
        ax10.step(np.arange(0,nf),np.nanmedian(T,axis=0),where='mid')        
        ax01.step(np.nanmedian(T,axis=1),np.arange(0,nsubint),where='mid')
        ax00.tick_params(labelbottom=False)
        ax01.tick_params(labelleft=False)
        ax00.set_ylabel('Subintegration')
        ax01.set_xlabel('OFF')
        ax10.set_xlabel('Channel')#,fontweight='bold')
        ax10.set_ylabel('OFF')
        ax10.set_ylim(10,40)
        ax00.set_title(calfits,fontsize=5)
        cax = fig.add_axes([ax00.get_position().x0,ax00.get_position().y1,ax00.get_position().width,0.02])
        cb  = fig.colorbar(mappable=Tmap,cax=cax,pad=0,orientation='horizontal')
        cb.ax.tick_params(labelsize=15,direction='in',labeltop=True,labelbottom=False)
        cax.set_title(r'OFF',fontsize=15)
        plt.savefig('OFF2D_p{}.pdf'.format(i),format='pdf',bbox_inches='tight')
        plt.close()

    return  