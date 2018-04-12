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

def run_one(pid, file_name, imsz, wcs, csv_list, return_dict):
  count = np.zeros(imsz)
  step = 10
  tranges = []
  trange = [0, 0]

  hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%file_name)
  offsets = np.load('../data/%s/cata/offsets1_10_new.npy'%file_name)

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
  #trange = [0, 0]
  return_dict[pid] = (count, tranges)

if __name__ == '__main__':
  if False:
    skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    hdulist = pyfits.open('../data/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data

    step = 10
    offsets = np.load('../data/new/offsets300-1343_10_new.npy')
    length = offsets.shape[0]
    print length
    for i in range(length-1):
      offsets[i+1] = offsets[i+1] + offsets[i]
    offsets = np.append(offsets, np.array([[0.,0.]]), axis=0)
    print offsets.shape
    initial = int((co_data[300][0]-0.5)*1000)

    with open('../data/csv_list') as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          with open('../data/NUVphotons_5_cal_large_10_full.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile)
            data = []
            i = 1
            trange[0] = float(reader.next()[0])/1000
            id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
            for row in reader:
              time = int(row[0])
              if time/1000.>trange[0]+1342:
                print 'break'
                break
              row = cal_photon(offsets, initial, time, row, step)
              data.append(row[id])
              writer.writerow(row)
              if i%200000 == 0:
                H = hist(data, wcs, imsz)
                count += H
                trange[1] = float(data[-1][0])/1000
                data=H = None
                data = []
                gc.collect()
                print i
              i+=1
            if len(data)>0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_5_cal_large_10_sep_full.fits', clobber=False)


  if False:
    skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    skyrange = [24,15]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data

    step = 10
    offsets = np.load('../data/05_r/cata/offsets300-1343_10_new.npy')
    print offsets.shape
    initial = int((co_data[300][0]-0.5)*1000)

    with open('../data/csv_list05_r') as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          with open('../data/05_r/cata/NUVphotons_05_cal_cata_10_full.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile)
            data = []
            i = 1
            trange[0] = float(reader.next()[0])/1000
            id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
            for row in reader:
              time = int(row[0])
              if time/1000.>trange[0]+1342:
                print 'break'
                break
              row = cal_photon_r(offsets, initial, time, row, step)
              data.append(row[id])
              writer.writerow(row)
              if i%200000 == 0:
                H = hist(data, wcs, imsz)
                count += H
                trange[1] = float(data[-1][0])/1000
                data=H = None
                data = []
                gc.collect()
                print i
              i+=1
            if len(data)>0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_gPr_cata_10_sep_full_flip.fits', clobber=False)

#no cal
  if False:
    skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    skyrange = [24,15]
    #skypos, skyrange = ([0, 0], [1.33333336,1.33333336])
    tranges = []
    trange = [0, 0]
    #imsz = imagetools.deg2pix(skypos, skyrange, 0.0016666667)
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    #wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0016666667)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data

    initial = int((co_data[300][0]-0.5)*1000)

    with open('../data/csv_list') as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file) 
          data = []
          i = 1
          trange[0] = float(reader.next()[0])/1000
          id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
          for row in reader:
            time = int(row[0])
            if time/1000.>trange[0]+780:
              print 'break'
              break
            data.append(np.array(row,dtype='float64')[id])
            if i%200000 == 0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
              data=H = None
              data = []
              gc.collect()
              print i
            i+=1
          if len(data)>0:
            H = hist(data, wcs, imsz)
            count += H
            trange[1] = float(data[-1][0])/1000
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_scan.fits', clobber=False)
    #hdulist.writeto('../fits/count_map_05_gPr_sep_big.fits', clobber=False)

#focal plane 
  if False:
    skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    skyrange = [24,15]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data

    step = 10
    offsets = np.load('../data/05_r/cata/offsets300-1343_10_new.npy')
    print offsets.shape
    initial = int((co_data[300][0]-0.5)*1000)

    skypos_f, skyrange_f = ([0, 0], [1.33333336,1.33333336])
    imsz_f = imagetools.deg2pix(skypos_f, skyrange_f, 0.0016666667)
    count_f = np.zeros(imsz_f)
    wcs_f = imagetools.define_wcs(skypos_f,skyrange_f,width=False,height=False,verbose=0,pixsz=0.0016666667)
    print imsz_f

    cata = spi.load_obj('stars05')
    cata_a = np.array(list(cata))
    cata_len = len(cata)
    cata_a = cata_a*np.pi/180.
    X, Y, Z = pycoo.spherical_to_cartesian(1, cata_a[:,1], cata_a[:,0])
    cata_car = np.array([X,Y,Z], dtype='float64').T
    rad = np.cos(0.00250*np.pi/180.)

    with open('../data/csv_list05_full') as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          data = []
          data_f = []
          i = 1
          trange[0] = float(reader.next()[0])/1000
          id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
          id_f = np.array([1,0,0,0,0,0,1,1,0,0,0])>0
          for row in reader:
            time = int(row[0])
            if time/1000.>trange[0]+1342:
              print 'break'
              break
            row = cal_photon_r(offsets, initial, time, row, step)
            if time<initial:
              pass
            else:
              center = pycoo.spherical_to_cartesian(1, row[9]*np.pi/180., row[8]*np.pi/180.)
              sep = np.dot(cata_car, center)
              coo = cata_a[sep>=rad,:]
              if len(coo) == 0:
                data_f.append(row[id_f])
            data.append(row[id])
            if i%200000 == 0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
              data=H = None
              data = []

              print 'select:%d'%len(data_f)
              if len(data_f) != 0:
                H_f = hist(data_f, wcs_f, imsz_f)
                count_f += H_f
                data_f=H_f=None
              data_f=[]
              gc.collect()
              print i
            i+=1
          if len(data)>0:
            H = hist(data, wcs, imsz)
            count += H
            trange[1] = float(data[-1][0])/1000

            H_f = hist(data_f, wcs_f, imsz_f)
            count_f += H_f
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_gPr_cata_10_sep_full_latest_new1342.fits', clobber=False)

    hdu = pyfits.PrimaryHDU(count_f)
    hdu = imagetools.fits_header('NUV', skypos_f, tranges, skyrange_f, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_gPr_full_cata_device_r0_latest_nostar_new1342.fits', clobber=False)

#focal plane revised
  if False:
    scan_num = '05'
    scan_num = sys.argv[1]
    print scan_num
    tranges = []
    trange = [0, 0]
    skypos_f, skyrange_f = ([0, 0], [1.33333336,1.33333336])
    imsz_f = imagetools.deg2pix(skypos_f, skyrange_f, 0.0016666667)
    count_f = np.zeros(imsz_f)
    wcs_f = imagetools.define_wcs(skypos_f,skyrange_f,width=False,height=False,verbose=0,pixsz=0.0016666667)
    print imsz_f

    step = 10
    hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_000%s_0001-asprta.fits'%scan_num)
    co_data = hdulist[1].data
    offsets = np.load('../data/%s_new/cata/offsets1_10_new.npy'%scan_num)
    initial = int((co_data[1][0]-0.5)*1000)

    cata = spi.load_obj('star_set%s'%scan_num)
    cata_a = np.array(list(cata))
    cata_len = len(cata)
    cata_a = cata_a*np.pi/180.
    X, Y, Z = pycoo.spherical_to_cartesian(1, cata_a[:,1], cata_a[:,0])
    cata_car = np.array([X,Y,Z], dtype='float64').T
    rad = np.cos(0.00250*np.pi/180.)

    with open('../data/csv_list%s_full'%scan_num) as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          data = []
          data_f = []
          i = 1
          trange[0] = float(reader.next()[0])/1000
          id_f = np.array([1,0,0,0,0,0,1,1,0,0,0])>0
          for row in reader:
            time = int(row[0])
            row = cal_photon_r(offsets, initial, time, row, step)
            if False:#time<initial:
              pass
            else:
              center = pycoo.spherical_to_cartesian(1, row[9]*np.pi/180., row[8]*np.pi/180.)
              sep = np.dot(cata_car, center)
              coo = cata_a[sep>=rad,:]
              if len(coo) == 0:
                data_f.append(row[id_f])
            if i%200000 == 0:
              print 'select:%d'%len(data_f)
              if len(data_f) != 0:
                H_f = hist(data_f, wcs_f, imsz_f)
                count_f += H_f
                data_f=H_f=None
              data_f=[]
              gc.collect()
              print i
            i+=1
          if len(data)>0:
            H = hist(data, wcs, imsz)
            count += H
            trange[1] = float(data_f[-1][0])/1000

            H_f = hist(data_f, wcs_f, imsz_f)
            count_f += H_f
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    hdu = pyfits.PrimaryHDU(count_f)
    hdu = imagetools.fits_header('NUV', skypos_f, tranges, skyrange_f, hdu=hdu, wcs=wcs_f)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_%s_gPr_cata_10_device_r0_nostar.fits'%scan_num, clobber=False)


  if False:
    #skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    skypos = [267.5, -26.]
    skyrange = [25,20]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    step = 10
    scan_list = ['05','14','23','32','41','50','59','68']

    with open('../data/csv_list05_14_new') as f:
      csv_list = f.read().splitlines() 
      file_num = 0
      print file_num
      for csv_file in csv_list:
        scan_num = scan_list[file_num]
        print scan_num
        hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_000%s_0001-asprta.fits'%scan_num)
        offsets = np.load('../data/%s_new/cata/offsets1_10_new.npy'%scan_num)
        '''
        if file_num == 0:
          hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
          offsets = np.load('../data/05_new/cata/offsets300-1343_10_new.npy')
          time_len = 1342
          print offsets.shape
        elif file_num == 1:
          hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00014_0001-asprta.fits')
          offsets = np.load('../data/14_new/cata/offsets300-1352_10_new.npy')
          time_len = 1351
          print offsets.shape
        '''
        co_data = hdulist[1].data
        initial = int((co_data[1][0]-0.5)*1000)

        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          data = []
          i = 1
          trange[0] = float(reader.next()[0])/1000
          id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
          for row in reader:
            time = int(row[0])
            '''
            if time/1000.>trange[0]+time_len:
              print 'break'
              break
            '''
            row = cal_photon_r(offsets, initial, time, row, step)
            data.append(row[id])
            if i%200000 == 0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
              data=H = None
              data = []
              gc.collect()
              print i
            i+=1
          if len(data)>0:
            H = hist(data, wcs, imsz)
            count += H
            trange[1] = float(data[-1][0])/1000
          file_num += 1
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05-68_gPr_cata_10_cut.fits', clobber=False)


#process
  if True:
    #logging.basicConfig(filename='../output/example.log',level=logging.DEBUG)
    #logging.info('start')

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

      asprta = np.load('../data/photon_list/%s_asp.npy'%name)
      sky_data = SkyCoord(asprta[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
      gal = sky_data.transform_to(Galactic)
      asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
      gal_l.append(np.mean(asprta[:,0]))

    skypos = [np.mean(gal_l), 0.]
    skyrange = [np.max(gal_l)-np.min(gal_l)+2., 24.]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix_g(skypos, skyrange, 0.002)
    wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    count = np.zeros(imsz)
    print skypos, skyrange
    print imsz

    step = 10

    total_p_num = 10
    manager = Manager()
    return_dict = manager.dict()
    p_list = []

    pid = 0
    for name in name_list:
      csvs = csv_lists[name]
      num_csv = len(csvs)
      p_num = int(round(1.*num_csv/total_num_csv*total_p_num))
      p_num = p_num if p_num>0 else 1
      print 'p_num:%d'%p_num
      chunk_len = int(num_csv/p_num)
      chunks=[csvs[x:x+100] for x in xrange(0, num_csv, chunk_len)]

      for chunk in chunks:
        p = Process(target=run_one, args=(pid, name, imsz, wcs, chunk, return_dict))
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
    hdulist.writeto('../fits/count_map_%s_gPr_cata_10_gal.fits'%name_file, clobber=False)

