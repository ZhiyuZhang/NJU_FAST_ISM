import numpy as np
from read_cal import read_cal
from read_spec import read_spec
from astropy.time import Time as ATime
from astropy.io import fits
import sys

sys.path.append('/Leo_disk/lingrui/my_projects/FAST_explore/FISH-pipeline/FISH_scripts')

# --------------- extract cal
# ls -d -1 $PWD/*fits > mylist
from extract_crafts_cal import extract_crafts_cal

psrlist = np.loadtxt('mylist',dtype=str)
for i in psrlist:
    extract_crafts_cal(i)

# ---------------- generate the gain cal table
from gaincal import gaincal
gaincal(callist='mycallist',outname='mycal',plot=True)




# ------------------ cal backend
calpath   = 'HighCAL/'
calfits1  = 'Dec+2150-M01_0001_210204_135403.fits'
calhdu1   = read_cal(calpath+calfits1)
calfits2  = 'Dec+2150-M01_0002_210204_135415.fits'
calhdu2   = read_cal(calpath+calfits2)
calfits3  = 'Dec+2150-M01_0003_210204_135428.fits'
calhdu3   = read_cal(calpath+calfits2)

calmjd1 = calhdu1[2].data['mjd']
calut1  = ATime(calmjd1, format='mjd').datetime

# ------------------ spec backend
specpath  = 'HIdata/'
specfits1 = 'Dec+2150_06_05_arcdrift-M01_N_0001.fits'
spechdu1  = read_spec(specpath+specfits1)
specfits2 = 'Dec+2150_06_05_arcdrift-M01_N_0002.fits'
spechdu2  = read_spec(specpath+specfits2)

specmjd1 = spechdu1[1].data['utobs'] 
specut1  = ATime(specmjd1, format='mjd').datetime

