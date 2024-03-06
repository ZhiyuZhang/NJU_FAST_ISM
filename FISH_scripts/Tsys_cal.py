# Generate the Tsys_psr

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

def timebin(onoff,cal_mjd,n):
    onoff_rebin   = []
    cal_mjd_rebin = []
    nsubint = len(cal_mjd)
    nbin    = nsubint//n
    for i in range(nbin):
        #print(i*n,i*n+n)
        onoff_rebin.append(np.nanmedian(onoff[i*n:i*n+n],axis=0))
        cal_mjd_rebin.append(np.nanmedian(cal_mjd[i*n:i*n+n]))
        # maybe it doesn't matter whether we choose mean/median in time?
    if nsubint%n != 0:
        onoff_rebin.append(np.nanmedian(onoff[-(nsubint%n):],axis=0))
        cal_mjd_rebin.append(np.nanmedian(cal_mjd[-(nsubint%n):]))
    onoff_rebin   = np.array(onoff_rebin)
    cal_mjd_rebin = np.array(cal_mjd_rebin)
    return onoff_rebin,cal_mjd_rebin

def gaincal(callist,outname='Tsys_psr',do_cal=True,Tcalfits='Tcal/20201014/CAL.20201014.low.W.fits',beam=1,sm_nf=16,sm_nt=300,plot=False):
    
    callist = np.loadtxt(callist,dtype=str)
    
    cal_mjd = []
    cal_on  = []
    cal_off = []
    
    for calfile in callist:
        try:
            hdu0,hdu1,hdu2 = read_crafts_cal(calfile)
            cal_mjd.append(hdu2.data['mjd'])
            cal_on.append(hdu2.data['cal'][:,1,:,:])
            cal_off.append(hdu2.data['cal'][:,0,:,:])
        except:
            print('WARNING: BAD FITS FILE!! Please check '+ calfile)
    
    cal_mjd = np.concatenate(cal_mjd,axis=0)
    cal_on  = np.concatenate(cal_on,axis=0)
    cal_off = np.concatenate(cal_off,axis=0)
    onoff   = (cal_on-cal_off)/20480*(20480+28672)/2
    nsubint = len(cal_mjd)
    
    del cal_on

    obsfreq  = hdu0.header['OBSFREQ']                # observing frequency
    obsnchan = hdu0.header['OBSNCHAN']               # observing channel number
    obsbw    = hdu0.header['OBSBW']                  # observing band width
    fmin     = float(obsfreq - obsbw/2.)             # minimum frequency
    fmax     = float(obsfreq + obsbw/2.)             # maximum frequency
    nf       = obsnchan                              # channel number
    df       = hdu1.header['CHAN_BW']                # channel width
    cal_freq = np.arange(fmin,fmax,df)/1e3           # frequency array; MHz to GHz
    p        = onoff.shape[1]                        # n_polarization
    
    onoff_interp      = interp(cal_mjd    , onoff      , cal_mjd  , kind = 'linear' , axis = 0 , smooth = sm_nt , smoothmode = 'nearest') # smooth along mjd
    
    if do_cal:
        
        _,noisehdu1 = read_crafts_spec(Tcalfits)
        noise_data  = np.transpose(noisehdu1.data['tcal'][0,beam-1,:,:],axes=(1,0))
        print('There are {} polarizations in the Tcal noise.'.format(noise_data.shape[0]))
        if noise_data.shape[0]<p: 
        # add NaN array when the length of the polarization axis of the noise data is smaller than 4.
            noise_data = np.concatenate((noise_data,np.nan*np.ones(shape=(p-noise_data.shape[0],noise_data.shape[1]))),axis=0)
        noise_freq  = noisehdu1.data['freq'][0]/1e3 # MHz to GHz

        noise_data_interp = interp(noise_freq , noise_data , cal_freq , kind = 'linear' , axis = 1 , smooth = sm_nf  , smoothmode = 'nearest') # smooth along frequency

        Tpsr              = noise_data_interp/onoff_interp*cal_off ### Tsys_psr
        #onoff_rebin,cal_mjd_rebin=timebin(onoff,cal_mjd,n)

        TSYS_PSR = fits.Column(name='TSYS_PSR', format="%dE" % (p*nf), array=Tpsr,dim=str((nf, p)))
        TIME     = fits.Column(name='MJD', format="D", array=cal_mjd)
        cols     = fits.ColDefs([TSYS_PSR,TIME])
        tsystab  = fits.BinTableHDU.from_columns(cols)
        hduout   = fits.HDUList([hdu0,hdu1,tsystab])

        cal_date = ATime(cal_mjd,format='mjd').datetime
        newname = outname + '-M{}'.format(str(beam).rjust(2,'0')) + '_' + cal_date[0].strftime("%y%m%d") + "_" + cal_date[0].strftime("%H%M%S") + 'to' \
                  + cal_date[-1].strftime("%y%m%d") + "_" + cal_date[-1].strftime("%H%M%S")
        hduout.writeto(newname+'.fits', overwrite=True)
        hduout.close()

    # ------------------------- plot on-off and resampled noise spectra
    if plot:
        for i in range(2):
            NoiInj        = onoff[:,i,:]
            NoiInj_interp = onoff_interp[:,i,:]
            fig = plt.figure()
            gs  = fig.add_gridspec(nrows=2,ncols=2,hspace=0,wspace=0,height_ratios=np.array([2,1]),width_ratios=np.array([2,1]))
            ax00   = fig.add_subplot(gs[0,0])
            ax01   = fig.add_subplot(gs[0,1],sharey=ax00)
            ax10   = fig.add_subplot(gs[1,0],sharex=ax00)
            Noimap = ax00.imshow(NoiInj,cmap='jet',aspect='auto',vmin=0,vmax=1)
            ax10.step(np.arange(0,nf),np.nanmedian(NoiInj,axis=0),where='mid')
            ax01.step(np.nanmedian(NoiInj,axis=1),np.arange(0,nsubint),where='mid')
            ax10.step(np.arange(0,nf),np.nanmedian(NoiInj_interp,axis=0),where='mid')
            ax01.step(np.nanmedian(NoiInj_interp,axis=1),np.arange(0,nsubint),where='mid')
            ax00.tick_params(labelbottom=False)
            ax01.tick_params(labelleft=False)
            ax00.set_ylabel('Subintegration')
            ax01.set_xlabel('ON-OFF')
            ax10.set_xlabel('Channel')#,fontweight='bold')
            ax10.set_ylabel('ON-OFF')
            cax = fig.add_axes([ax00.get_position().x0,ax00.get_position().y1,ax00.get_position().width,0.02])
            cb  = fig.colorbar(mappable=Noimap,cax=cax,pad=0,orientation='horizontal')
            cb.ax.tick_params(labelsize=15,direction='in',labeltop=True,labelbottom=False)
            cax.set_title(r'$Readout$',fontsize=15)
            plt.savefig('onoff2D_p{}.pdf'.format(i),format='pdf',bbox_inches='tight')
            plt.close()


        if do_cal:
            
            fig  = plt.figure()
            gs   = fig.add_gridspec(nrows=2,ncols=1,hspace=0,height_ratios=np.array([1,1]))
            ax00 = fig.add_subplot(gs[0,0])
            ax10 = fig.add_subplot(gs[1,0],sharex=ax00)
            ax00.step(noise_freq,noise_data[0],where='mid',label='AA')
            ax00.step(cal_freq,noise_data_interp[0],linestyle='--',where='mid',label='AA resampled')
            ax10.step(noise_freq,noise_data[1],where='mid',label='BB')
            ax10.step(cal_freq,noise_data_interp[1],linestyle='--',where='mid',label='BB resampled')
            ax00.set_ylim(1.0,1.3)
            ax10.set_ylim(1.0,1.3)
            #ax00.set_xlim(1.40,1.44)
            ax00.legend()
            ax10.legend()
            ax10.set_xlabel('FREQ [GHz]')
            ax10.set_ylabel('T_CAL [K]')
            plt.savefig('Tcal_resample.pdf',format='pdf',bbox_inches='tight')
            plt.close()

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
                ax10.set_ylim(23,40)
                cax = fig.add_axes([ax00.get_position().x0,ax00.get_position().y1,ax00.get_position().width,0.02])
                cb  = fig.colorbar(mappable=Tmap,cax=cax,pad=0,orientation='horizontal')
                cb.ax.tick_params(labelsize=15,direction='in',labeltop=True,labelbottom=False)
                cax.set_title(r'$T_{\rm A}$ [K]',fontsize=15)
                plt.savefig('Tpsr2D_p{}.pdf'.format(i),format='pdf',bbox_inches='tight')
                plt.close()

    return  

if __name__ == '__main__':
    callist = sys.argv[1]
    plot = sys.argv[2]
    gaincal(callist,plot)
