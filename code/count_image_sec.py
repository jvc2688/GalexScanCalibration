import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc
import h5py
from astropy.io import fits as pyfits
import pos_range

if __name__ == '__main__':
  initial_sec = 300
  hdulist = pyfits.open('../data/AIS_GAL_SCAN_00005_0001-asprta.fits')
  co_data = hdulist[1].data
  intitial_asp = co_data[initial_sec]
  skypos = [intitial_asp[1], intitial_asp[2]]
  skyrange = [1.38,1.38]

  tranges = []
  trange = [0, 0]
  imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
  count = np.zeros(imsz)
  print imsz
  print skypos
  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)

  with open('../data/csv_list') as f:
    csv_list = f.read().splitlines() 
    for csv_file in csv_list:
      with open(csv_file, 'rb') as file:
        reader = csv.reader(file)
        data = []
        i = 1
        initial_time = trange[0] = int(intitial_asp[0]*1000)
        print initial_time
        for row in reader:
          time = int(row[0])
          if time < (initial_time + 0):
            continue
          if time >= (initial_time + 1000):
            print 'break'
            break
          else:
            data.append(row)
            trange[1] = time
            if i%100000 == 0:
              coo = np.array(data, dtype='float64')[:,-3:-1]
              foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
              H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                      bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
              count += H
              data=coo=foc=H = None
              gc.collect()
              data = []
            i+=1
        if len(data)>0:
          coo = np.array(data, dtype='float64')[:,-3:-1]
          foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
          H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                     bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
          count += H
        print i
      trange[0] = trange[0]/1000.
      trange[1] = trange[1]/1000.
      tranges.append(trange)
      trange = [0, 0]
  print tranges
  print count.shape
  hdu = pyfits.PrimaryHDU(count)
  hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
  hdulist = pyfits.HDUList([hdu])
  hdulist.writeto('../fits/count_map_5_sec%d_%d_try.fits'%(initial_sec, initial_sec+1), clobber=False)

