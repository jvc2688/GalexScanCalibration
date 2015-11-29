import matplotlib
matplotlib.use('Agg')
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

def run_cal(name):
  asp_file = np.load('../data/photon_list/%s_asp.npy'%name)
  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
  offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%name)

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

  np.save('../data/photon_list/%s_asp_cal.npy'%name, asp_file)

def interpolate_offsets(name, interval=0.5, offsets=None):
  dt = 0.005
  begin = interval/2/dt
  middle = interval/dt-1
  end = interval/2/dt-1

  asp_file = np.load('../data/photon_list/%s_asp.npy'%name)
  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
  if offsets==None:
    offsets = np.load('../data/%s/cata/offsets1_10_new_half.npy'%name)
  print middle
  mask = np.ones(offsets.shape[0])
  print offsets.shape[0]*100
  insert = np.zeros((middle,2))
  insert_mask = np.zeros(middle)
  for i in range(1, offsets.shape[0]):
    mask = np.insert(mask, middle*(i-1)+i, insert_mask, axis=0)
    offsets = np.insert(offsets, middle*(i-1)+i, insert, axis=0)
  insert = np.zeros((begin,2))
  insert_mask = np.zeros(begin)
  mask = np.insert(mask, 0, insert_mask, axis=0)
  offsets = np.insert(offsets, 0, insert, axis=0)
  insert = np.zeros((end,2))
  insert_mask = np.zeros(end)
  offsets = np.concatenate([offsets, insert], axis=0)
  mask = np.concatenate([mask, insert_mask], axis=0)
  mask = mask>0
  index = np.arange(offsets.shape[0])
  for i in range(2):
    offsets[~mask, i] = np.interp(index[~mask], index[mask], offsets[mask,i])

  print offsets.shape

  np.save('../data/%s/cata/offsets1_10_new_inter_half_fine.npy'%name, offsets)

  f, axes = plt.subplots(2, 1, squeeze=False)
  axes[0,0].plot(offsets[:,0], '.b')
  axes[0,0].set_ylabel('RA')
  axes[1,0].plot(offsets[:,1], '.b')
  axes[1,0].set_ylabel('DEC')
  plt.setp(axes[0,0].get_xticklabels(), visible=False)
  plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                   wspace=0, hspace=0)
  plt.savefig('../plots/%s/cata/offsets1_10_new_inter_half_fine.png'%name, dpi=190)
  plt.clf()

  co_data = hdulist[1].data
  initial = int((co_data[1][0]-0.5)*1000)

  step = 200

  t0 = asp_file[0,0]
  for i in range(asp_file.shape[0]):
    row = asp_file[i,:]
    time = int(row[0]*1000)
    asp_file[i,:] = cal_photon_r(offsets, initial, time, row, step)

  np.save('../data/photon_list/%s_asp_cal_inter_half.npy'%name, asp_file)

def secondary_cal(name):
  asp_file = np.load('../data/photon_list/%s_asp.npy'%name)
  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
  offsets = np.load('../data/%s/cata/offsets1_10_new_inter_half_fine.npy'%name)
  offsets_sec = np.load('../data/%s/cata/offsets1_10_new_sec.npy'%name)

  for i in range(offsets_sec.shape[0]):
    offsets[i*200:(i+1)*200] += offsets_sec[i]

  np.save('../data/%s/cata/offsets_inter_half_sec.npy'%name, offsets)

  co_data = hdulist[1].data
  initial = int((co_data[1][0]-0.5)*1000)

  step = 200

  t0 = asp_file[0,0]
  for i in range(asp_file.shape[0]):
    row = asp_file[i,:]
    time = int(row[0]*1000)
    asp_file[i,:] = cal_photon_r(offsets, initial, time, row, step)

  np.save('../data/photon_list/%s_asp_cal_inter_half_sec.npy'%name, asp_file)
  print 'asp caled'

if __name__ == '__main__':
  
  if False:
    name = 'AIS_GAL_SCAN_00257_0001'
    interpolate_offsets(name)


  #origin
  if False:
    name = sys.argv[1]
    asp_file = np.load('../data/photon_list/%s_asp.npy'%name)
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%name)

    co_data = hdulist[1].data
    initial = int((co_data[1][0]-0.5)*1000)

    step = 10

    t0 = asp_file[0,0]
    for i in range(asp_file.shape[0]):
      row = asp_file[i,:]
      time = int(row[0]*1000)

      asp_file[i,:] = cal_photon_r(offsets, initial, time, row, step)

    np.save('../data/photon_list/%s_asp_cal.npy'%name, asp_file)

#interpolate
  if False:
    interval = 0.5
    begin = interval/2/0.005
    middle = interval/0.005-1
    end = interval/2/0.005-1
    name = sys.argv[1]
    asp_file = np.load('../data/photon_list/%s_asp.npy'%name)
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    offsets = np.load('../data/%s/cata/offsets1_10_new_half.npy'%name)
    mask = np.ones(offsets.shape[0])
    print offsets.shape[0]*100
    insert = np.zeros((middle,2))
    insert_mask = np.zeros(middle)
    for i in range(1, offsets.shape[0]):
      print i
      mask = np.insert(mask, middle*(i-1)+i, insert_mask, axis=0)
      offsets = np.insert(offsets, middle*(i-1)+i, insert, axis=0)
    insert = np.zeros((begin,2))
    insert_mask = np.zeros(begin)
    mask = np.insert(mask, 0, insert_mask, axis=0)
    offsets = np.insert(offsets, 0, insert, axis=0)
    insert = np.zeros((end,2))
    insert_mask = np.zeros(end)
    offsets = np.concatenate([offsets, insert], axis=0)
    mask = np.concatenate([mask, insert_mask], axis=0)
    mask = mask>0
    print mask[0:60]
    index = np.arange(offsets.shape[0])
    for i in range(2):
      offsets[~mask, i] = np.interp(index[~mask], index[mask], offsets[mask,i])

    print offsets.shape

    np.save('../data/%s/cata/offsets1_10_new_inter_half_fine.npy'%name, offsets)

    co_data = hdulist[1].data
    initial = int((co_data[1][0]-0.5)*1000)

    plt.plot(offsets[:,1],'.b')
    plt.show()

    step = 200

    t0 = asp_file[0,0]
    for i in range(asp_file.shape[0]):
      row = asp_file[i,:]
      time = int(row[0]*1000)

      #if time/1000. > t0+1351:
      #  print 'break'
      #  break

      asp_file[i,:] = cal_photon_r(offsets, initial, time, row, step)

    np.save('../data/photon_list/%s_asp_cal_inter_half.npy'%name, asp_file)

#interpolate
  if True:
    interval = 0.1
    begin = interval/2/0.005
    middle = interval/0.005-1
    end = interval/2/0.005-1
    name = sys.argv[1]
    asp_file = np.load('../data/photon_list/%s_asp.npy'%name)
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%name)
    mask = np.ones(offsets.shape[0])
    print offsets.shape[0]*20
    insert = np.zeros((middle,2))
    insert_mask = np.zeros(middle)
    for i in range(1, offsets.shape[0]):
      #print i
      mask = np.insert(mask, middle*(i-1)+i, insert_mask, axis=0)
      offsets = np.insert(offsets, middle*(i-1)+i, insert, axis=0)
    insert = np.zeros((begin,2))
    insert_mask = np.zeros(begin)
    mask = np.insert(mask, 0, insert_mask, axis=0)
    offsets = np.insert(offsets, 0, insert, axis=0)
    insert = np.zeros((end,2))
    insert_mask = np.zeros(end)
    offsets = np.concatenate([offsets, insert], axis=0)
    mask = np.concatenate([mask, insert_mask], axis=0)
    mask = mask>0
    print mask[0:60]
    index = np.arange(offsets.shape[0])
    for i in range(2):
      offsets[~mask, i] = np.interp(index[~mask], index[mask], offsets[mask,i])

    print offsets.shape

    np.save('../data/%s/cata/offsets_inter_half_sec.npy'%name, offsets)

    co_data = hdulist[1].data
    initial = int((co_data[1][0]-0.5)*1000)

    plt.plot(offsets[:,1],'.b')
    plt.show()

    step = 200

    t0 = asp_file[0,0]
    for i in range(asp_file.shape[0]):
      row = asp_file[i,:]
      time = int(row[0]*1000)

      #if time/1000. > t0+1351:
      #  print 'break'
      #  break

      asp_file[i,:] = cal_photon_r(offsets, initial, time, row, step)

    np.save('../data/photon_list/%s_asp_cal_inter_half_sec.npy'%name, asp_file)

  if False: 
    name = sys.argv[1]
    asp_file = np.load('../data/photon_list/%s_asp.npy'%name)
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    offsets = np.load('../data/%s/cata/offsets1_10_new_inter_half_fine.npy'%name)
    offsets_sec = np.load('../data/%s/cata/offsets1_10_new_sec.npy'%name)

    for i in range(offsets_sec.shape[0]):
      offsets[i*200:(i+1)*200] += offsets_sec[i]

    np.save('../data/%s/cata/offsets_inter_half_sec.npy'%name, offsets)

    co_data = hdulist[1].data
    initial = int((co_data[1][0]-0.5)*1000)

    step = 200

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

    np.save('../data/photon_list/%s_asp_cal_inter_half_sec.npy'%name, asp_file)


