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
from scipy import interpolate
import gnomonic as gn
import h5py

def dis_correct(data, dis_map, asp):
    ix = np.digitize(data['time']/1000., asp[:,0])-1

    low_mask = data['xa']<16
    high_mask = data['xa']>=16

    coo = np.array([data['xi'][low_mask],data['eta'][low_mask]]).T
    data['xi'][low_mask] -= dis_map['xi_l'](coo)
    data['eta'][low_mask] -= dis_map['eta_l'](coo)

    coo = np.array([data['xi'][high_mask],data['eta'][high_mask]]).T
    data['xi'][high_mask] -= dis_map['xi_h'](coo)
    data['eta'][high_mask] -= dis_map['eta_h'](coo)

    ra, dec = gn.gnomrev_simple(data['xi'], data['eta'], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return np.array([ra,dec]).T

def detector2gon(pos_list, size):
  pos_array = np.array(pos_list, dtype='float64')
  pos = (pos_array/50.-0.5)*36000*800*0.001666
  return pos

def get_key(path):
  name = os.path.basename(path)
  return int(re.split('\.', name)[0])

def split_seq(seq, size):
  newseq = []
  splitsize = 1.0/size*len(seq)
  for i in range(size):
          newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
  return newseq

def distr_p_num(csv_info_list, total_p_num, total_num_csv):
  total_p_num = total_p_num - csv_info_list.shape[0]
  p_used = 0
  for csv_info in csv_info_list:
    num = csv_info['num']
    p_num = int(1.*num/total_num_csv*total_p_num)
    csv_info['p_num'] += p_num
    p_used += p_num

  p_remain = total_p_num-p_used
  if p_remain > 0:
    csv_info_list = np.sort(csv_info_list, order='num')[::-1]
    for i in range(p_remain):
      csv_info_list[i]['p_num'] += 1
  return csv_info_list


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

def cal_photon_fast(offsets, initial, data, step):
  off_len = offsets.shape[0]
  num = data.shape[0]
  index_list = np.zeros(num)
  corr = np.zeros((num,2))
  index = np.floor((data[:,0]-initial)*step/1000.).astype(int)
  mask = np.logical_and(index<off_len, data[:,0]>=initial)
  index = index[mask]
  corr[mask,:] = offsets[index]
  data[:,1:3] -= corr
  return data[:,1:3]

def hist(data, wcs, imsz):
  coo = np.array(data, dtype='float64')[:,1:]
  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

def hist_g(data, wcs, imsz):
  #coo = np.array(data, dtype='float64')#[:,1:]
  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(data,1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H


def run_one_inter_fast_f5(pid, file_name, imsz, wcs, csv_list, suffix, asp, scst_time, hv, dis_map, return_dict):

  count = np.zeros(imsz)
  step = 200
  tranges = []
  trange = [0, 0]

  #hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/{0}{1}-asprta.fits'.format(file_name, suffix))
  '''
  scst_file = pyfits.open('../AIS_GAL_SCAN/scst/{0}-scst.fits'.format(file_name))
  scst_data = scst_file[1].data
  scst_time = scst_data['pktime'].copy()
  hv = scst_data['hvnom_nuv'].copy()
  scst_file.close()
  '''
  f = h5py.File('../data/{0}-cal-sec/photon_list.hdf5'.format(file_name),'r')
  col = ('time', 'xa', 'q', 'xi', 'eta')
  data = f['data'][col]
  f.close()

  print '{0}:loaded,shape:{1}'.format(pid,data.shape)
  #trim the data with scst flag
  scst_ix = np.digitize(data['time']/1000., scst_time)
  ix_mask = scst_ix<scst_time.shape[0]
  scst_ix = scst_ix[ix_mask]
  data = data[ix_mask]
  data = data[hv[scst_ix]>0]
  print '{0}:filter,shape:{1}'.format(pid,data.shape)
  if len(data)==0:
    print '{0}:skip'.format(pid)

  trange[0] = float(data['time'][0])/1000.
  trange[1] = float(data['time'][-1])/1000.

  data = dis_correct(data, dis_map, asp)
  print '{0}:correct,shape:{1}'.format(pid,data.shape)
  print '{0}:convert2gal'.format(pid)
  sky_data = SkyCoord(data, unit='deg', frame=FK5, equinox='J2000.0')
  gal = sky_data.transform_to(Galactic)
  data = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
  print '{0}:hist'.format(pid) 
  H = hist_g(data, wcs, imsz)
  count += H
  sky_data=gal=data=H = None
  gc.collect()
  print '{0}:{1}'.format(pid,i)
  tranges.append(trange)
  #return_dict[pid] = (count, tranges)
  np.save('../fits/scan_map/%s_gal_sec_count_tmp%d.npy'%(name_file, pid), count)
  return_dict[pid] = tranges

if __name__ == '__main__':

#process interpolate fast
  if True:
    #logging.basicConfig(filename='../output/example.log',level=logging.DEBUG)
    #logging.info('start')
    pixsz = 0.000416666666666667 #0.00166667
    name_file = sys.argv[1]
    if len(sys.argv)>3:
      suffix = sys.argv[3]
    else:
      suffix = ''

    with open('../name_scan/%s'%name_file) as f:
      name_list = f.read().splitlines()
    print name_list

    csv_lists = {}
    csv_info_list = []
    total_num_csv = 0
    gal_l = []
    dtype1 = np.dtype([('name', np.str_, 23), ('num', int), ('p_num', int)])
    for name in name_list:
      csvs = glob.glob("../data/{0}{1}/split/*.csv".format(name, suffix))
      csv_lists[name] = csvs
      num = len(csvs)
      csv_info_list.append((name, num, 1))
      total_num_csv += num

      asprta = np.load('../data/photon_list/{0}{1}_asp.npy'.format(name, suffix))#np.load('../data/photon_list/%s_asp.npy'%name)
      sky_data = SkyCoord(asprta[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
      gal = sky_data.transform_to(Galactic)
      asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
      print np.min(asprta[:,0])
      print np.max(asprta[:,0])
      print asprta
      gal_l.append(np.mean(asprta[:,0]))

    skypos = [np.mean(gal_l), 0.]
    print skypos
    print np.max(gal_l)
    print np.min(gal_l)
    skyrange = [np.max(gal_l)-np.min(gal_l)+2., 22.]
    print skyrange
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz)
    wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
    count = np.zeros(imsz)
    print skypos, skyrange
    print imsz

    detector_size = [50,50]

    dis_map_l = np.load('../data/distortion/5-185-xa_low-centroids.npy')
    dis_map_h = np.load('../data/distortion/5-185-xa_high-centroids.npy')

    dis_map_l = dis_map_l.reshape(50*50,2)/detector_size[0]*36000*800*0.001666
    dis_map_h = dis_map_h.reshape(50*50,2)/detector_size[0]*36000*800*0.001666
    cen_x = np.arange(detector_size[0])+0.5
    cen_y = np.arange(detector_size[1])+0.5
    xi = np.repeat(detector2gon(cen_x, detector_size[0]), detector_size[1])
    eta = np.tile(detector2gon(cen_y, detector_size[1]), detector_size[0])
    points = np.array([xi,eta]).T

    fxi_l = interpolate.LinearNDInterpolator(points, dis_map_l[:,0])
    feta_l = interpolate.LinearNDInterpolator(points, dis_map_l[:,1])
    fxi_h = interpolate.LinearNDInterpolator(points, dis_map_h[:,0])
    feta_h = interpolate.LinearNDInterpolator(points, dis_map_h[:,1])

    dis_map = {'xi_l':fxi_l, 'xi_h':fxi_h, 'eta_l':feta_l, 'eta_h':feta_h}

    csv_info_list = np.array(csv_info_list, dtype=dtype1)
    total_p_num = int(sys.argv[2])

    print csv_info_list
    print total_p_num
    print total_num_csv
    csv_info_list = distr_p_num(csv_info_list, total_p_num, total_num_csv)

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0

    print 'new'
    for csv_info in csv_info_list:
      name = csv_info['name']
      p_num = csv_info['p_num']
      num = csv_info['num']
      print num
      print 'p_num:%d'%p_num

      try:
        asp = np.load('../data/photon_list/{0}{1}_asp_new.npy'.format(name, suffix))
      except IOError:
        co_data = pyfits.open('../AIS_GAL_SCAN/asprta/{0}{1}-asprta.fits'.format(name, suffix))[1].data
        T = co_data['T']
        ra = co_data['ra']
        dec = co_data['dec']
        roll = co_data['roll']
        t_new = np.arange((T.shape[0]-1)*200)*0.005+T[0]
        ra_new = np.interp(t_new, T, ra)
        dec_new = np.interp(t_new, T, dec)
        roll_new = np.interp(t_new, T, roll)
        asp_solution_tmp = np.array([t_new, ra_new, dec_new, roll_new]).T
        np.save('../data/photon_list/{0}{1}_asp_new.npy'.format(name, suffix), asp)

      scst_file = pyfits.open('../AIS_GAL_SCAN/scst/{0}-scst.fits'.format(name))
      scst_data = scst_file[1].data
      scst_time = scst_data['pktime'].copy()
      hv = scst_data['hvnom_nuv'].copy()
      scst_file.close()

      csvs = csv_lists[name]
      csvs = sorted(csvs, key=get_key)

      chunks = split_seq(csvs, p_num)
      #print chunks
      for chunk in chunks:
        print len(chunk)
        p = Process(target=run_one_inter_fast_f5, args=(pid, name, imsz, wcs, chunk, suffix, asp, scst_time, hv, dis_map, return_dict))
        p.start()
        p_list.append(p)
        pid += 1
    print pid

    for p in p_list:
      p.join()
    print 'all done'

    for i in range(0, pid):
      count += np.load('../fits/scan_map/%s_gal_sec_count_tmp%d.npy'%(name_file, i)) #return_dict[i][0]
      os.remove('../fits/scan_map/%s_gal_sec_count_tmp%d.npy'%(name_file, i))
      tranges += return_dict[i]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/scan_map/count_map_{0}{1}_count_dis.fits'.format(name_file, suffix), clobber=False)

