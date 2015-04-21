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
  skypos = [268.73293, -29.573102]
  skyrange = [2., 2.]

  skypos, skyrange = pos_range.get_pos_range(name_list='../galex/name_list')
  tranges = []
  trange = [0, 0]
  imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
  count = np.zeros(imsz)
  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
  print imsz
  with open('../galex/csv_list') as f:
    csv_list = f.read().splitlines() 
    for csv_file in csv_list:
      with open(csv_file, 'rb') as file:
        reader = csv.reader(file)
        data = []
        i = 1
        trange[0] = float(reader.next()[0])/1000
        for row in reader:
          data.append(row)
          if i%100000 == 0:
            coo = np.array(data, dtype='float64')[:,-3:-1]
            foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
            H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                    bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
            count += H
            trange[1] = float(data[-1][0])/1000
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
  tbl = imagetools.movie_tbl('NUV', tranges, framesz=0)
  hdu = pyfits.PrimaryHDU(count)
  hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
  hdulist = pyfits.HDUList([hdu, tbl])
  hdulist.writeto('count_map_5_test1.fits', clobber=False)

