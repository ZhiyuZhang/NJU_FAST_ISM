#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits
import numpy as np
import os
import sys

dayfilename=sys.argv[1]
outnum=int(sys.argv[2])
# In[2]:


hdu_class=fits.open('oneobs.fits')
hdu_class[0].header
IC443_RA=94.51125
IC443_Dec=22.66


# In[ ]:


if not os.path.exists('divided_fits'):
    os.makedirs('divided_fits')
def makefits(hdu_fast,scan_num,ibeam):
    first_spec=np.zeros([1,1,1,hdu_fast[1].data['NCHAN'][scan_num]],dtype=np.float32)
    first_spec[0,0,0,:]=hdu_fast[1].data['data'][scan_num][:,0]

    newheader=hdu_class[0].header
    newheader['NAXIS1']=hdu_fast[1].data['NCHAN'][scan_num]
    newheader['DATAMIN']=np.min(first_spec)
    newheader['DATAMAX']=np.max(first_spec)
    newheader['BUNIT']='elec'
    newheader['CDELT1']=hdu_fast[1].data['CHAN_BW'][scan_num]*1e6
    newheader['CRPIX1']=0
    newheader['CRVAL1']=0
    newheader['CRVAL2']=hdu_fast[1].data['OBJ_RA'][scan_num]
    newheader['CRVAL3']=hdu_fast[1].data['OBJ_DEC'][scan_num]
    newheader['CDELT2']=hdu_fast[1].data['OBJ_RA'][scan_num]-IC443_RA
    newheader['CDELT3']=hdu_fast[1].data['OBJ_DEC'][scan_num]-IC443_Dec
    newheader['TELESCOP']='FAST'
    newheader['OBJECT']='IC443'
    newheader['LINE']='HI'
    newheader['RESTFREQ']=hdu_fast[1].data['FREQ'][scan_num]*1e6
    #newheader['DELTAV']=hdu_fast[1].data['CHAN_BW'][scan_num]*1e6/1.4204058e9*3e8
    #newheader['IMAGFREQ']=1.4204058e9
    newheader['TSYS']=1
    newheader['OBSTIME']=hdu_fast[1].data['EXPOSURE'][scan_num]
    newheader['SCAN-NUM']=scan_num+1
    newheader['DATE-OBS']=hdu_fast[1].data['DATE-OBS'][scan_num][0:-1]
    #newheader.remove('DATE-RED')
    if 'DELTAV' in newheader.keys():
        newheader.remove('DELTAV')
    if 'IMAGFREQ' in newheader.keys():
        newheader.remove('IMAGFREQ')
    if 'BLANK' in newheader.keys():
        newheader.remove('BLANK')

    hduout=fits.PrimaryHDU(first_spec,newheader)
    fitsname='./divided_fits/beam_%d_scan_%d.fits'%(ibeam,scan_num)
    os.system('rm '+fitsname)
    hduout.writeto(fitsname)
for i_obs in range(1,11):
    for i_beam in range(1,20):
        hdu_fast=fits.open('added_radec/'+dayfilename+'_arcdrift-M'+str(i_beam).rjust(2,'0')+'_N_00'+str(i_obs).rjust(2,'0')+'.fits')
        for i_scan in range(len(hdu_fast[1].data)):
            makefits(hdu_fast,i_scan,i_beam)
    os.system('rm class_fmt/obs_%d.fast'%(outnum+i_obs))
    os.system('class @load_fits.class %d'%(outnum+i_obs))


# In[ ]:




