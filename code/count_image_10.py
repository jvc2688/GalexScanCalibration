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
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
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
  off_len = offsets.shape[0]
  if time<initial:
    return np.array(row, dtype='float64')
  else:
    index = int(math.floor((time-initial)*step/1000.))
    row = np.array(row, dtype='float64')
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

def hist_g(data, wcs, imsz):
  coo = np.array(data, dtype='float64')#[:,1:]
  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

#interpolate and identify outliers
def generate_count(pid, file_name, imsz, wcs, csv_list, offsets_name, step, return_dict):
  count = np.zeros(imsz)
  tranges = []
  trange = [0, 0]

  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%file_name)
  offsets = np.load('../data/%s/cata/%s'%(file_name, offsets_name))

  co_data = hdulist[1].data
  initial = int((co_data[1][0]-0.5)*1000)

  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)
      data = []
      i = 1
      trange[0] = float(reader.next()[0])/1000
      id = np.array([0,0,0,0,0,0,0,0,1,1,0])>0
      for row in reader:
        time = int(row[0])
        row = cal_photon_r(offsets, initial, time, row, step)
        data.append(row[id])
        if i%200000 == 0:
          sky_data = SkyCoord(data, unit='deg', frame=FK5, equinox='J2000.0')
          gal = sky_data.transform_to(Galactic)
          data = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1) 
          H = hist_g(data, wcs, imsz)
          count += H
          trange[1] = float(time)/1000
          data=H = None
          data = []
          gc.collect()
          print i
        i+=1
      if len(data)>0:
        sky_data = SkyCoord(data, unit='deg', frame=FK5, equinox='J2000.0')
        gal = sky_data.transform_to(Galactic)
        data = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1) 
        H = hist_g(data, wcs, imsz)
        count += H
        trange[1] = float(time)/1000
  tranges.append(trange)
  return_dict[pid] = (count, tranges)

def generate_count_nocal(pid, file_name, imsz, wcs, csv_list, return_dict):
  count = np.zeros(imsz)
  step = 10
  tranges = []
  trange = [0, 0]

  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%file_name)

  co_data = hdulist[1].data
  initial = int((co_data[1][0]-0.5)*1000)

  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)
      data = []
      i = 1
      trange[0] = float(reader.next()[0])/1000
      id = np.array([0,0,0,0,0,0,0,0,1,1,0])>0
      for row in reader:
        time = int(row[0])
        row = np.array(row, dtype='float64')
        data.append(row[id])
        if i%200000 == 0:
          sky_data = SkyCoord(data, unit='deg', frame=FK5, equinox='J2000.0')
          gal = sky_data.transform_to(Galactic)
          data = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1) 
          H = hist_g(data, wcs, imsz)
          count += H
          trange[1] = float(time)/1000
          data=H = None
          data = []
          gc.collect()
          print i
        i+=1
      if len(data)>0:
        sky_data = SkyCoord(data, unit='deg', frame=FK5, equinox='J2000.0')
        gal = sky_data.transform_to(Galactic)
        data = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1) 
        H = hist_g(data, wcs, imsz)
        count += H
        trange[1] = float(time)/1000
  tranges.append(trange)
  return_dict[pid] = (count, tranges)


if __name__ == '__main__':
#generate count map, multiprocess
  if True:
    name_file = sys.argv[1]
    with open('../name_new/%s'%name_file) as f:
      name_list = f.read().splitlines()
    print name_list

    csv_lists = {}
    total_num_csv = 0
    gal_l = []
    for name in name_list:
      csvs = glob.glob("../data/%s/split/*.csv"%name)
      csv_lists[name] = csvs
      total_num_csv += len(csvs)

      asprta = np.load('../data/photon_list/%s_asp_cal_inter_half_sec.npy'%name)#np.load('../data/photon_list/%s_asp.npy'%name)
      sky_data = SkyCoord(asprta[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
      gal = sky_data.transform_to(Galactic)
      asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
      gal_l.append(np.mean(asprta[:,0]))

    skypos = [np.mean(gal_l), 0.]
    skyrange = [np.max(gal_l)-np.min(gal_l)+2., 24.]
    tranges = []
    trange = [0, 0]
    pixsz = 00166667
    imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz)
    wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
    count = np.zeros(imsz)
    print skypos, skyrange
    print imsz

    total_p_num = int(sys.argv[2])
    manager = Manager()
    return_dict = manager.dict()
    p_list = []

    step = 200
    offsets_name = 'offsets_inter_half_sec.npy'

    pid = 0
    for name in name_list:
      csvs = csv_lists[name]
      csvs = sorted(csvs, key=get_key)
      print len(csvs)
      num_csv = len(csvs)
      p_num = int(round(1.*num_csv/total_num_csv*total_p_num))
      p_num = p_num if p_num>0 else 1
      print 'p_num:%d'%p_num

      chunks = split_seq(csvs, p_num)
      for chunk in chunks:
        print len(chunk)
        p = Process(target=generate_count, args=(pid, name, imsz, wcs, chunk, offsets_name, step, return_dict))
        p.start()
        p_list.append(p)
        pid += 1
    print pid

    for p in p_list:
      p.join()
    print 'all done'

    for i in range(0, pid):
      count += return_dict[i][0]
      tranges += return_dict[i][1]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_%s_gPr_cata_10_gal_inter_half_sec.fits'%(name_file), clobber=False)

 #process, no CAL
  if False:
    name_file = sys.argv[1]
    with open('../name_new/%s'%name_file) as f:
      name_list = f.read().splitlines()
    print name_list

    csv_lists = {}
    total_num_csv = 0
    gal_l = []
    for name in name_list:
      csvs = glob.glob("../data/%s/split/*.csv"%name)
      csv_lists[name] = csvs
      total_num_csv += len(csvs)

      asprta = np.load('../data/photon_list/%s_asp.npy'%name)#np.load('../data/photon_list/%s_asp.npy'%name)
      sky_data = SkyCoord(asprta[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
      gal = sky_data.transform_to(Galactic)
      asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
      gal_l.append(np.mean(asprta[:,0]))

    pixsz = 0.002
    skypos = [np.mean(gal_l), 0.]
    skyrange = [np.max(gal_l)-np.min(gal_l)+2., 24.]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz)
    wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
    count = np.zeros(imsz)
    print skypos, skyrange
    print imsz

    step = 10

    total_p_num = int(sys.argv[2])
    manager = Manager()
    return_dict = manager.dict()
    p_list = []

    pid = 0
    for name in name_list:
      csvs = csv_lists[name]
      csvs = sorted(csvs, key=get_key)
      csvs = csvs[sec:sec+1]
      print len(csvs)
      num_csv = len(csvs)
      p_num = int(round(1.*num_csv/total_num_csv*total_p_num))
      p_num = p_num if p_num>0 else 1
      print 'p_num:%d'%p_num

      chunks = split_seq(csvs, p_num)
      #print chunks
      for chunk in chunks:
        print len(chunk)
        p = Process(target=run_one_nocal, args=(pid, name, imsz, wcs, chunk, return_dict))
        p.start()
        p_list.append(p)
        pid += 1
    print pid

    for p in p_list:
      p.join()
    print 'all done'

    for i in range(0, pid):
      count += return_dict[i][0]
      tranges += return_dict[i][1]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_%s_gPr_oncal_gal.fits'%(name_file), clobber=False)

