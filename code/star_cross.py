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
import os
import glob
from multiprocessing import Process
from multiprocessing import Manager
import logging
import re

def get_key(path):
  name = os.path.basename(path)
  return int(re.split('\.', name)[0])

def split_seq(seq, size):
  newseq = []
  splitsize = 1.0/size*len(seq)
  for i in range(size):
    newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
  return newseq

def cal_photon(offsets, initial, time, row, step):
  if time<initial:
    return np.array(row, dtype='float64')
  else:
    index = int(math.floor((time-initial)*step/1000.-1))
    row = np.array(row, dtype='float64')
    row[-3:-1] += offsets[index]
    return row 

def cal_photon_r(offsets, initial, time, row, step):
  row = np.array(row, dtype='float64')
  if row.shape[0]<11:
    print 'error: %d'%time
  off_len = offsets.shape[0]
  if time<initial:
    return row
  else:
    index = int(math.floor((time-initial)*step/1000.))
    if index >= off_len:
      return row
    row[-3:-1] -= offsets[index]
    return row 

def hist(data, wcs, imsz):
  coo = np.array(data, dtype='float64')[:,1:]
  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

def hist_flat(data, imsz):
  coo = np.array(data, dtype='float64')[:,1:]
  foc = ((coo/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  #H,xedges,yedges=np.histogram2d(foc[:,0], foc[:,1],\
                            #bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

def run_one(pid, file_name, imsz_list, wcs_list, cata_list, csv_list, return_dict):
  count_list = []
  for imsz in imsz_list:
    count = np.zeros(imsz)
    count_list.append(count)

  star_num = len(count_list)

  cata_list[:,1:3] = cata_list[:,1:3]*np.pi/180.
  X, Y, Z = pycoo.spherical_to_cartesian(1, cata_list[:,2], cata_list[:,1])
  cata_car = np.array([X,Y,Z], dtype='float64').T
  aperture_size = 8./60./60.
  rad = np.cos(aperture_size*np.pi/180.)
  annulus = np.cos([9./60./60.*np.pi/180., 11./60./60.*np.pi/180.])

  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%file_name)
  offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%file_name)

  co_data = hdulist[1].data
  initial = int((co_data[1][0]-0.5)*1000)

  step = 10

  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)
      data = []
      i = 1
      trange[0] = float(reader.next()[0])/1000
      id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
      for row in reader:
        time = int(row[0])
        row = cal_photon_r(offsets, initial, time, row, step)
        #add row into data
        #data.append(row[id])
        center = pycoo.spherical_to_cartesian(1, row[9]*np.pi/180., row[8]*np.pi/180.)
        sep = np.dot(cata_car, center)
        coo = cata_list[sep>=rad,:]
        ann_coo = cata_list[np.logical_and(sep>=annulus[1], sep<=annulus[0]),:]
        for photon in coo:
          with open('../data/%s/star/extra/star_%d-%d.csv'%(file_name, photon[0], pid), 'ab') as star_csv:
            writer = csv.writer(star_csv)
            writer.writerow(row)
          star_csv.close()

        for photon in ann_coo:
          with open('../data/%s/star/extra/star_bkg_%d-%d.csv'%(file_name, photon[0], pid), 'ab') as star_bkg_csv:
            writer = csv.writer(star_bkg_csv)
            writer.writerow(row)
          star_bkg_csv.close()

        '''
        #hist the data
        if i%200000 == 0:
          for star in range(star_num):
            H = hist(data, wcs_list[star], imsz_list[star])
            count_list[star] += H
          trange[1] = float(time)/1000
          data=H = None
          data = []
          gc.collect()
          print i
        i+=1
      #hist the remaining data
      if len(data)>0:
        for star in range(star_num):
          H = hist(data, wcs_list[star], imsz_list[star])
          count_list[star] += H
        trange[1] = float(time)/1000
      '''
  tranges.append(trange)
  return_dict[pid] = (count_list, tranges)

def run_one_flat(pid, file_name, imsz, cata_list, csv_list, return_dict):
  count = np.zeros(imsz)

  #star_num = len(count_list)

  cata_list[:,1:3] = cata_list[:,1:3]*np.pi/180.
  X, Y, Z = pycoo.spherical_to_cartesian(1, cata_list[:,2], cata_list[:,1])
  cata_car = np.array([X,Y,Z], dtype='float64').T
  aperture_size = 8./60./60.
  rad = np.cos(aperture_size*np.pi/180.)
  annulus = np.cos([9./60./60.*np.pi/180., 11./60./60.*np.pi/180.])

  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%file_name)
  offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%file_name)

  co_data = hdulist[1].data
  initial = int((co_data[1][0]-0.5)*1000)

  step = 10

  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)
      data = []
      i = 1
      trange[0] = float(reader.next()[0])/1000
      id = np.array([1,0,0,0,0,0,1,1,0,0,0])>0
      for row in reader:
        time = int(row[0])
        row = cal_photon_r(offsets, initial, time, row, step)
        #add row into data
        #data.append(row[id])
        center = pycoo.spherical_to_cartesian(1, row[9]*np.pi/180., row[8]*np.pi/180.)
        sep = np.dot(cata_car, center)
        coo = cata_list[sep>=rad,:]
        #ann_coo = cata_list[np.logical_and(sep>=annulus[1], sep<=annulus[0]),:]
        if len(coo)>0:
          data.append(row[id])
        '''
        for photon in coo:
          with open('../data/%s/star/extra/star_%d-%d.csv'%(file_name, photon[0], pid), 'ab') as star_csv:
            writer = csv.writer(star_csv)
            writer.writerow(row)
          star_csv.close()
        '''
        '''
        for photon in ann_coo:
          with open('../data/%s/star/extra/star_bkg_%d-%d.csv'%(file_name, photon[0], pid), 'ab') as star_bkg_csv:
            writer = csv.writer(star_bkg_csv)
            writer.writerow(row)
          star_bkg_csv.close()
        '''

        #hist the data
        if i%200000 == 0:
          H = hist_flat(data, imsz)
          count += H
          data=H = None
          data = []
          gc.collect()
          print i
        i+=1
      #hist the remaining data
      if len(data)>0:
        H = hist_flat(data, imsz)
        count += H

  return_dict[pid] = count

def run_flat_chunks(pid, file_name, imsz, cata_list, csv_list, return_dict):
  count = np.zeros(imsz)

  #star_num = len(count_list)
  trange = [0,0]

  cata_list[:,1:3] = cata_list[:,1:3]*np.pi/180.
  X, Y, Z = pycoo.spherical_to_cartesian(1, cata_list[:,2], cata_list[:,1])
  cata_car = np.array([X,Y,Z], dtype='float64').T
  aperture_size = 8./60./60.
  rad = np.cos(aperture_size*np.pi/180.)
  annulus = np.cos([9./60./60.*np.pi/180., 11./60./60.*np.pi/180.])

  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%file_name)
  offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%file_name)

  co_data = hdulist[1].data
  initial = int((co_data[1][0]-0.5)*1000)

  step = 10
  with open(csv_list[0], 'rb') as file:
    reader = csv.reader(file)
    trange[0] = int(reader.next()[0])

  time = 0
  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)
      data = []
      i = 1
      id = np.array([1,0,0,0,0,0,1,1,0,0,0])>0
      for row in reader:
        time = int(row[0])

        row = cal_photon_r(offsets, initial, time, row, step)
        #add row into data
        #data.append(row[id])
        center = pycoo.spherical_to_cartesian(1, row[9]*np.pi/180., row[8]*np.pi/180.)
        sep = np.dot(cata_car, center)
        coo = cata_list[sep>=rad,:]
        #ann_coo = cata_list[np.logical_and(sep>=annulus[1], sep<=annulus[0]),:]
        if len(coo)>0:
          data.append(row[id])

        '''
        for photon in coo:
          with open('../data/%s/star/extra/star_%d-%d.csv'%(file_name, photon[0], pid), 'ab') as star_csv:
            writer = csv.writer(star_csv)
            writer.writerow(row)
          star_csv.close()
        '''
        '''
        for photon in ann_coo:
          with open('../data/%s/star/extra/star_bkg_%d-%d.csv'%(file_name, photon[0], pid), 'ab') as star_bkg_csv:
            writer = csv.writer(star_bkg_csv)
            writer.writerow(row)
          star_bkg_csv.close()
        '''

        #hist the data
        if i%200000 == 0:
          H = hist_flat(data, imsz)
          count += H
          data=H = None
          data = []
          gc.collect()
          print i
        i+=1
      #hist the remaining data
      if len(data)>0:
        H = hist_flat(data, imsz)
        count += H

  trange[1] = time
  return_dict[pid] = (count, trange)

if __name__ == '__main__':

  if False:
    name = sys.argv[1]
    #name = 'AIS_GAL_SCAN_00032_0001'
    print name

    output = '../data/%s/star/extra/star.csv'%(name)
    path = os.path.dirname(output)
    if not os.path.exists(path):
      os.makedirs(path)
    else:
      print 'exists'

    cata = spi.load_obj('../data/%s_starset_extra_full'%name)
    cata_a = np.array(list(cata))
    cata_len = len(cata)
    cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)

    count_list, imsz_list, wcs_list = [],[],[]
    skypos_list, skyrange_list = [],[]
    for i in range(cata_list.shape[0]):
      skypos, skyrange = ([cata_list[i,1], cata_list[i,2]] , [0.03, 0.03])
      print skypos
      print skyrange
      skypos_list.append(skypos)
      skyrange_list.append(skyrange)
      imsz = imagetools.deg2pix(skypos, skyrange, 0.000416666666666667)
      count = np.zeros(imsz)
      wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.000416666666666667)
      count_list.append(count)
      wcs_list.append(wcs)
      imsz_list.append(imsz)

    step = 10
    tranges = []
    trange = [0, 0]

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%name)
    initial = int((co_data[1][0]-0.5)*1000)

    csvs = glob.glob("../data/%s/split/*.csv"%name)

    total_p_num = 19
    manager = Manager()
    return_dict = manager.dict()
    p_list = []

    pid = 0
    csvs = sorted(csvs, key=get_key)
    chunks = split_seq(csvs, total_p_num)
    for chunk in chunks:
      print len(chunk)
      p = Process(target=run_one, args=(pid, name, imsz_list, wcs_list, cata_list, chunk, return_dict))
      p.start()
      p_list.append(p)
      pid += 1
    print pid

    for p in p_list:
      p.join()
    print 'all done'
    '''
    for i in range(0, pid):
      for j in range(0, cata_list.shape[0]):
        count_list[j] += return_dict[i][0][j]
      tranges += return_dict[i][1]

    print tranges
    for i in range(0, cata_list.shape[0]):
      hdu = pyfits.PrimaryHDU(count_list[i])
      hdu = imagetools.fits_header('NUV', skypos_list[i], tranges, skyrange_list[i], hdu=hdu, wcs=wcs_list[i])
      hdulist = pyfits.HDUList([hdu])
      hdulist.writeto('../fits/%s/extra/count_map_star%d.fits'%(name,i), clobber=False)
    '''

#flat field---star photons
  if False:
    name = sys.argv[1]
    #name = 'AIS_GAL_SCAN_00014_0001'
    print name

    output = '../fits/%s/extra/flat.fits'%(name)
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
      os.makedirs(dir)
    else:
      print 'exists'

    cata = spi.load_obj('../data/%s_starset_extra_all_star'%name)
    cata_a = np.array(list(cata))
    cata_a = cata_a[np.logical_and(cata_a[:,4]>15.5,cata_a[:,4]<17.5)]
    cata_len = len(cata_a)
    print cata_len
    cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)

    imsz = [800, 800]
    count = np.zeros(imsz)
    pixsz = 0.0016666667

    step = 10
    tranges = []
    trange = [0, 0]

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%name)
    initial = int((co_data[1][0]-0.5)*1000)

    csvs = glob.glob("../data/%s/split/*.csv"%name)

    total_p_num = 19
    manager = Manager()
    return_dict = manager.dict()
    p_list = []

    pid = 0
    csvs = sorted(csvs, key=get_key)
    chunks = split_seq(csvs, total_p_num)
    for chunk in chunks:
      print len(chunk)
      p = Process(target=run_one_flat, args=(pid, name, imsz, cata_list, chunk, return_dict))
      p.start()
      p_list.append(p)
      pid += 1
    print pid

    for p in p_list:
      p.join()
    print 'all done'
  
    for i in range(0, pid):
      count += return_dict[i]

    flat = pyfits.open('../data/cal/NUV_flat.fits')

    hdu = pyfits.PrimaryHDU(count)
    hdu.header = flat[0].header
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/%s/extra/flat.fits'%(name), clobber=False)

#flat field---star photons, chunks
  if True:
    name = sys.argv[1]
    #name = 'AIS_GAL_SCAN_00014_0001'
    print name

    output = '../fits/%s/extra/new/flat.fits'%(name)
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
      os.makedirs(dir)
    else:
      print 'exists'

    cata = spi.load_obj('../data/%s_starset_extra_all_star'%name)
    cata_a = np.array(list(cata))
    cata_a = cata_a[np.logical_and(cata_a[:,4]>15.5,cata_a[:,4]<17.5)]
    cata_len = len(cata_a)
    print cata_len
    cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)

    imsz = [800, 800]
    count = np.zeros(imsz)
    pixsz = 0.0016666667

    step = 10
    tranges = []

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%name)
    initial = int((co_data[1][0]-0.5)*1000)

    csvs = glob.glob("../data/%s/split/*.csv"%name)

    total_p_num = 20
    manager = Manager()
    return_dict = manager.dict()
    p_list = []

    pid = 0
    csvs = sorted(csvs, key=get_key)
    chunks = split_seq(csvs, total_p_num)
    for chunk in chunks:
      print len(chunk)
      p = Process(target=run_flat_chunks, args=(pid, name, imsz, cata_list, chunk, return_dict))
      p.start()
      p_list.append(p)
      pid += 1
    print pid

    for p in p_list:
      p.join()
    print 'all done'
    
    flat = pyfits.open('../data/cal/NUV_flat.fits')
    tranges = []
    for i in range(0, pid):
      count, trange = return_dict[i]

      hdu = pyfits.PrimaryHDU(count)
      hdu.header = flat[0].header
      hdulist = pyfits.HDUList([hdu])
      hdulist.writeto('../fits/%s/extra/new/flat%d.fits'%(name, i), clobber=False)
      tranges.append(trange)

    np.save('../fits/%s/extra/new/tranges.npy'%name, tranges)

