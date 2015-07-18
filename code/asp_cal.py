import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc
from astropy.io import fits as pyfits
import pos_range
import math
import spilt_csv as spi
import astropy.coordinates as pycoo
import sys

def cal_photon_r(offsets, initial, time, row, step):
  off_len = offsets.shape[0]
  if time<initial:
    return row
  else:
    index = int(math.floor((time-initial)*step/1000.))
    if index >= off_len:
      return row
    row[1:3] -= offsets[index]
    return row 

if __name__ == '__main__':

  scan_num = sys.argv[1]

  asp_file = np.load('../data/photon_list/NUVPhoton%s_full_new_asp.npy'%scan_num)
  hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_000%s_0001-asprta.fits'%scan_num)
  offsets = np.load('../data/%s_new/cata/offsets1_10_new.npy'%scan_num)

  co_data = hdulist[1].data
  initial = int((co_data[1][0]-0.5)*1000)

  step = 10

  t0 = asp_file[0,0]
  for i in range(asp_file.shape[0]):
    row = asp_file[i,:]
    time = int(row[0]*1000)
    '''
    if time/1000. > t0+1351:
      print 'break'
      break
    '''
    asp_file[i,:] = cal_photon_r(offsets, initial, time, row, step)

  np.save('../data/photon_list/NUVPhoton%s_full_new_asp_cal_cata_wa.npy'%scan_num, asp_file)
