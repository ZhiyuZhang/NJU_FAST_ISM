import sys
import xlrd
import re
import os
import numpy as np
from datetime import datetime,timedelta
from astropy.table import Table
from astropy.io import fits
import ephem
import time

tablename=sys.argv[1]
dir_name=sys.argv[2]

xyztable = xlrd.open_workbook(tablename)

#compute center ra dec from x,y,z

convert=np.pi/180.0

class coordinate:
      def __init__(place):
          place.name=''
          place.ra=''
          place.dec=''
          place.alt=''

Dawodang      = coordinate()
Dawodang.name = 'Dawodang'
Dawodang.ra   = str(106.85645657571428)
Dawodang.dec  = str(25.6534387)
Dawodang.alt  = str(1138.72178)

def localtoec(az,el,obstime):
    getlocal                    = ephem.Observer()
    getlocal.long, getlocal.lat = Dawodang.ra,Dawodang.dec # at Dawodang
    getlocal.pressure           = 0
    getlocal.epoch              = ephem.J2000
    CurrentUT                   = time.localtime(obstime)
    CurrentUT3                  = time.strftime("%Y/%m/%d %H:%M:%S",CurrentUT)
    getlocal.date               = CurrentUT3 
    ra,dec = getlocal.radec_of(az*convert,el*convert)
    return ra,dec

def xyztoradec(xyz,obstime):
    R=np.sqrt(xyz[0]**2+xyz[1]**2+xyz[2]**2)
    z=-xyz[2]
    y=-xyz[0]
    x=-xyz[1]
    az0 = 0.0
    if (np.abs(x) < 1.0e-8):                 # If X is a very small number
       if (y>0):
          az0 = 90.0                         # abs(X) <1E-8; Y >0;  AZ = 90
       else:
          az0 = 270.0                        # abs(X) <1E-8; Y <=0; AZ = 270
    else:
       tempaz0 = np.arctan(y/x)/convert      # X >=1E-8;  arctan(y/x) to degree
       if (x>0 and y>0):
          az0 = tempaz0                      # X > 0; Y > 0; AZ = tempaz0
       if ((x>0 and y<0) or (x>0 and y==0)):
           az0 = tempaz0+360.0               # X > 0; Y <= 0;   AZ = tempaz0+360
       if (x<0 and y>0):
           az0 = tempaz0+180.0               # X < 0; Y > 0;    AZ = tempaz0+180
       if ((x<0 and y<0) or (x<0 and y==0)):
           az0 = tempaz0+180.0               # X < 0; Y <= 0;   AZ = tempaz0+180
    el0 = np.arcsin(z/R)/convert             # elevation to degree
    ra,dec=localtoec(az0,el0,obstime)
    return ra/convert,dec/convert
ra,dec=xyztoradec([0.008134610019624,-29.6900997161865,-159.009002685547],1612077097.851)
ra,dec




#generate position calibration table
table = xyztable.sheet_by_name('整控-馈源舱数据')
result_table_210203=np.zeros([table.nrows-1,4])
for i in range(table.nrows-1):
    result_table_210203[i,0]=(datetime.fromisoformat(table.cell_value(i+1,0))+timedelta(hours=-8)).timestamp()
    result_table_210203[i,1:3]=xyztoradec(table.row_values(i+1)[7:10],result_table_210203[i,0])
    result_table_210203[i,3]=table.cell_value(i+1,18)



files=os.listdir(dir_name)
fitsfiles=[]
narrow_fits=[]
beamswiftorigin={'01':(0,0),'02':(5.74,0.00811),'03':(2.88,-4.97),'04':(-2.86,-4.98),'05':(-5.74,-0.0127),
                '06':(-2.88,4.97),'07':(2.86,4.98),'08':(11.5,0.0116),'09':(8.61,-4.96),'10':(5.75,-9.93),
                '11':(0.018,-9.94),'12':(-5.71,-9.96),'13':(-8.6,-4.99),'14':(-11.5,-0.03),'15':(-8.63,4.95),
                '16':(-5.77,9.93),'17':(-0.0181,9.94),'18':(5.73,9.95),'19':(8.61,4.98)}
for i in range(len(files)):
    app=re.search(r'.([a-z|A-Z]*?)$',files[i]).group(0)
    if app == '.fits':
        fitsfiles.append(files[i])
        if re.search('_N_',files[i])!=None:
            narrow_fits.append(files[i])



for fitsfile in narrow_fits:
    print('Transforming: '+fitsfile)
    dirfile=dir_name+'/'+fitsfile
    hdu=fits.open(dirfile)
    beam_num=re.findall('-M(\d+)_',fitsfile)[0]
    start_search_i=0
    for obs in range(len(hdu[1].data)):
        obstime= datetime.fromisoformat(hdu[1].data['DATE-OBS'][obs].replace('T',' ').replace('Z','')).timestamp()
        if obstime<result_table_210203[0,0] or obstime>result_table_210203[-1,0]:
            print('数据与表格不匹配！')
            exit()
        for i in range(start_search_i,len(result_table_210203)):    # To find where observation time in cal table
            if obstime>=result_table_210203[i,0] and obstime<result_table_210203[i+1,0]:
                time_pos_i=i
                start_search_i=i
                break
        left_value=result_table_210203[time_pos_i]
        right_value=result_table_210203[time_pos_i+1]
        true_value=left_value+(obstime-left_value[0])/(right_value[0]-left_value[0])*(right_value-left_value)
        sin_theta=true_value[3] # Here to control anti- or clockwise. Additional '-' means clockwise 
        cos_theta=np.sqrt(1-true_value[3]**2)
        rotat_matrix=np.array([[cos_theta,-sin_theta],[sin_theta,cos_theta]])
        beamswift=np.dot(rotat_matrix,np.array(beamswiftorigin[beam_num])) # Unit is arcmin
        RA=true_value[1]+beamswift[0]/60
        Dec=true_value[2]+beamswift[1]/60
        hdu[1].data['OBJ_RA'][obs]=RA
        hdu[1].data['OBJ_DEC'][obs]=Dec
    os.system('rm added_radec/'+fitsfile)
    hdu.writeto('added_radec/'+fitsfile)

