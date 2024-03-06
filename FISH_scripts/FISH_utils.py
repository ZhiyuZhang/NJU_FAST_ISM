import numpy as np
from scipy.interpolate import interp1d
from scipy import ndimage

def interp(x,y,xnew,kind='linear',axis=0,smooth=0,smoothmode='nearest'):
    if smooth != 0:
        y = ndimage.uniform_filter1d(y,smooth,axis=axis,mode=smoothmode)
    interp_func = interp1d(x=x,y=y,kind=kind,axis=axis,fill_value='extrapolate')
    ynew = interp_func(xnew)
    return ynew

def rebin(array,rebin_n=2):

    if rebin_n<=1:
        return array

    array_trim  = array[0:int(len(array)/rebin_n)*rebin_n] # drop a tail of the array such that it can be grouped every n iterms

    array_rebin = np.nanmean(np.array(array_trim).reshape(-1,rebin_n),axis=1)

    return array_rebin

