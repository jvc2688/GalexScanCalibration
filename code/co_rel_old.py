from __future__ import print_function
import sys
import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc
from astropy.io import fits as pyfits
import catalog_fits
import astropy.coordinates as pycoo


def get_corr_map(coo1, coo2, skypos, skyrange):
  imsz = imagetools.deg2pix(skypos, skyrange, 0.0001)
  count = np.zeros(imsz)
  print(imsz)
  co_rel = np.array([[0,0]])
  len1 = coo1.shape[0]
  len2 = coo2.shape[0]
  print(len1,len2)
  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)
  if len2>len1:
    for i in range(len2):
      print(i)
      co_rel = np.concatenate((co_rel, np.roll(coo2, i, axis=0)[0:len1,:]-coo1), axis = 0)
      if (i+1)%200 == 0:
        foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
        H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                   bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
        count += H
        co_rel = np.array([[0, 0]])
  else:
    for i in range(len1):
      print(i)
      co_rel = np.concatenate((co_rel, coo2-np.roll(coo1, i, axis=0)[0:len2,:]), axis = 0)
      if (i+1)%200 == 0:
        foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
        H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                   bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
        count += H
        co_rel = np.array([[0, 0]])
  if co_rel.shape[0]>1:
    foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
    H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                               bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
    count += H
  return count

def get_data(tranges):
  data = []
  trange = [0 , 0]
  num_interval = len(tranges)
  tranges = np.array(tranges, dtype='double')
  tranges = (tranges*1000).astype(int)
  with open('../data/csv_list') as f:
    csv_list = f.read().splitlines() 
    for csv_file in csv_list:
      with open(csv_file, 'rb') as file:
        reader = csv.reader(file)
        for i in range(num_interval):
          data.append([])
        initial_time = trange[0] = int(tranges[0,0]) #int(reader.next()[0]) + 0
        final_time = trange[1] = int(tranges[num_interval-1, 1])
        print(initial_time)
        for row in reader:
          time = int(float(row[0]))
          if time < (initial_time + 0):
            continue
          if time >= (final_time):
            print('break')
            break
          for i in range(num_interval):
            if time >= tranges[i, 0] and time < tranges[i, 1]:
              data[i].append(row)
              break
  return data

def angle_filter(data, skypos, angle):
  center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
  rad = np.cos(angle*np.pi/180.)
  print(center, rad)

  data_rad = data*np.pi/180.
  X, Y, Z = pycoo.spherical_to_cartesian(1, data_rad[:,1], data_rad[:,0])
  data_car = np.array([X,Y,Z], dtype='float64').T

  sep = np.dot(data_car, center)
  coo = data[sep>=rad, :]

  return coo

if __name__ == '__main__':
  initial_sec = int(sys.argv[1])

  hdulist = pyfits.open('../data/AIS_GAL_SCAN_00005_0001-asprta.fits')
  co_data = hdulist[1].data
  intitial_asp = co_data[initial_sec]
  center = np.array([intitial_asp[1], intitial_asp[2]])  

  print(intitial_asp)

  skypos = [0.0, 0.0]
  skyrange = [0.04, 0.04]

  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)


  initial_time = intitial_asp[0]
  tranges = [[initial_time, initial_time+1], [initial_time+1, initial_time+2]]
  print(tranges)

  data = get_data(tranges)

  coo1 = np.array(data[0], dtype='float64')[:,-3:-1]
  coo2 = np.array(data[1], dtype='float64')[:,-3:-1]

  coo1 = angle_filter(coo1, center, 0.4)
  coo2 = angle_filter(coo2, center, 0.4)

  #point_pos, time = catalog_fits.get_pos_time('../Galex/AIS_GAL_SCAN_00005_0001-asprta.fits', 452)
  #coo2 = catalog_fits.get_catalog(point_pos, 0.69)

  count = get_corr_map(coo2, coo1, skypos, skyrange)

  hdu = pyfits.PrimaryHDU(count)
  hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
  hdulist = pyfits.HDUList([hdu])
  hdulist.writeto('../fits/co/co_map%d_%d_zoom_nocal.fits'%(initial_sec, initial_sec+1), clobber=False)
    
