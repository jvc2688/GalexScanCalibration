import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc
from astropy.io import fits as pyfits
import pos_range
import math

def cal_photon(offsets, initial, time, row):
  if time<initial:
    return row
  else:
    index = int(math.floor((time-initial)/1000.-1))
    row = np.array(row, dtype='float64')
    row[-3:-1] += offsets[index]
    #print initial, time, index
    return row 

if __name__ == '__main__':
  skypos = [268.73293, -29.573102]
  skyrange = [2., 2.]

  skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
  tranges = []
  trange = [0, 0]
  imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
  count = np.zeros(imsz)
  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
  print imsz

  hdulist = pyfits.open('../data/AIS_GAL_SCAN_00005_0001-asprta.fits')
  co_data = hdulist[1].data

  offsets = np.load('../data/offsets300-1342_r.npy')
  length = offsets.shape[0]
  for i in range(length-1):
    offsets[i+1] = offsets[i+1] + offsets[i]
  np.append(offsets, np.array([0.,0.]))
  initial = int(co_data[300][0]*1000)

  with open('../data/csv_list') as f:
    csv_list = f.read().splitlines() 
    for csv_file in csv_list:
      with open(csv_file, 'rb') as file:
        reader = csv.reader(file)
        with open('../data/NUVphotons_5_cal_whole.csv', 'wb') as csvfile:
          writer = csv.writer(csvfile)
          data = []
          i = 1
          trange[0] = float(reader.next()[0])/1000
          for row in reader:
            time = int(row[0])
            row = cal_photon(offsets, initial, time, row)
            data.append(row)
            writer.writerow(row)
            if i%100000 == 0:
              coo = np.array(data, dtype='float64')[:,-3:-1]
              foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
              H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                      bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
              count += H
              trange[1] = float(data[-1][0])/1000
              if trange[1]>trange[0]+790:
                break
              data=coo=foc=H = None
              data = []
              gc.collect()
              print i
            i+=1
          if len(data)>0:
            coo = np.array(data, dtype='float64')[:,-3:-1]
            foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
            H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                       bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
            count += H
            trange[1] = float(data[-1][0])/1000
      tranges.append(trange)
      trange = [0, 0]

  print tranges
  print count.shape
  hdu = pyfits.PrimaryHDU(count)
  hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
  hdulist = pyfits.HDUList([hdu])
  hdulist.writeto('../fits/count_map_5_cal_new_sep.fits', clobber=False)

