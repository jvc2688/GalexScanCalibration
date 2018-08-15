import matplotlib
#matplotlib.use('Agg')
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
import centroid_kde as ck
import cPickle as pickle
import catalog_fits
import os
from multiprocessing import Process
from multiprocessing import Manager
import spilt_csv as spi
import math
import centroid_plot_csv_new as centroid_csv
import gnomonic as gn
import c3
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
import subprocess


def get_detector_pos(pos_list, size):
  pos_array = np.array(pos_list, dtype='float64')
  pos = ((pos_array/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*size
  return pos

def cart2pol(x, y):
  rho = np.sqrt(x**2 + y**2)
  phi = np.arctan2(y, x)
  return rho, phi

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def hist_g(data, wcs, imsz):
  coo = np.array(data, dtype='float64')#[:,1:]
  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

def get_corr_map(coo1, coo2, hist_size):
  co_rel = np.array([[0,0]])
  len1 = coo1.shape[0]
  len2 = coo2.shape[0]
  #print(len1,len2)
  if len2>len1:
    for i in range(len1):
      #print(i)
      tmp_co = coo2-coo1[i,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=0.2,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=0.2,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
  else:
    for i in range(len2):
      #print(i)
      tmp_co = coo2[i,:]-coo1
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=0.2,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=0.2,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)

  H,xedges,yedges=np.histogram2d(co_rel[1:,1], co_rel[1:,0],\
                             bins=hist_size, range=([ [-0.2,0.2],[-0.2,0.2] ]))

  return H


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

def cal_asp_r(offsets, initial, time, row, step):
  off_len = offsets.shape[0]
  if time<initial:
    return np.array(row, dtype='float64')
  else:
    index = int(math.floor((time-initial)*step/1000.))
    row = np.array(row, dtype='float64')
    if index >= off_len:
      return row
    row[0:2] -= offsets[index]
    return row 

def get_data(tranges, scan_name):
  look_up = load_obj('../data/%s-cal-sec/split/lookup'%scan_name)
  data = []
  trange = [0 , 0]
  num_interval = len(tranges)
  tranges = np.array(tranges, dtype='double')
  tranges = (tranges*1000).astype(int)

  initial_time = trange[0] = int(tranges[0,0]) #int(reader.next()[0]) + 0
  final_time = trange[1] = int(tranges[num_interval-1, 1])
  for i in range(200):
    if '%d'%(initial_time) in look_up:
      print 'initial found:%d'%initial_time
      break
    initial_time += 5

  for i in range(200):
    if '%d'%(final_time) in look_up:
      print 'final found:%d'%final_time
      break
    final_time -= 5

  for i in range(num_interval):
    data.append([])

  if initial_time>=final_time:
    return data
  
  start_index = look_up['%d'%initial_time]
  end_index = look_up['%d'%final_time]
  csv_list = []
  for i in range(start_index, end_index+1):
    csv_list.append('../data/%s-cal-sec/split/%d.csv'%(scan_name, i))
 
  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)

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

def get_data_star(tranges, scan_name):
  look_up = load_obj('../data/%s-cal-sec-star/split/lookup'%scan_name)
  data = []
  trange = [0 , 0]
  num_interval = len(tranges)
  tranges = np.array(tranges, dtype='double')
  tranges = (tranges*1000).astype(int)

  initial_time = trange[0] = int(tranges[0,0]) #int(reader.next()[0]) + 0
  final_time = trange[1] = int(tranges[num_interval-1, 1])
  for i in range(200):
    if '%d'%(initial_time) in look_up:
      print 'initial found:%d'%initial_time
      break
    initial_time += 5

  for i in range(200):
    if '%d'%(final_time) in look_up:
      print 'final found:%d'%final_time
      break
    final_time -= 5

  for i in range(num_interval):
    data.append([])

  if initial_time>=final_time:
    return data
  
  start_index = look_up['%d'%initial_time]
  end_index = look_up['%d'%final_time]
  csv_list = []
  for i in range(start_index, end_index+1):
    csv_list.append('../data/%s-cal-sec-star/split/%d.csv'%(scan_name, i))
 
  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)

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
  #print(center, rad)

  data_rad = data*np.pi/180.
  X, Y, Z = pycoo.spherical_to_cartesian(1, data_rad[:,1], data_rad[:,0])
  data_car = np.array([X,Y,Z], dtype='float64').T

  sep = np.dot(data_car, center)
  coo = data[sep>=rad, :]

  return coo

def angle_filter_out(data, skypos, angle):
  center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
  rad = np.cos(angle*np.pi/180.)
  #print(center, rad)

  data_rad = data*np.pi/180.
  X, Y, Z = pycoo.spherical_to_cartesian(1, data_rad[:,1], data_rad[:,0])
  data_car = np.array([X,Y,Z], dtype='float64').T

  sep = np.dot(data_car, center)
  coo = data[sep<=rad, :]

  return coo

def run_one_r_sec(pid, scan_name, step, asp_cal, start, end, hist_size, detector_size, nbin, return_dict):

    print('run one r sec')

    num_co = int(1/step)
    
    #load asp file
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]
    print length

    #hist_size = [50,50]
    #detector_size = [100,100]
    bins = []
    for bin_num in range(nbin*2):
      hist_list = []
      for i in range(detector_size[0]):
        for j in range(detector_size[1]):
          hist_list.append(np.zeros(hist_size))
      bins.append(hist_list)

    for initial_sec in range(start, end):
      print 'num:%d'%initial_sec

      intitial_asp = co_data[initial_sec]
      initial_time = intitial_asp[0]-step/2.

      tranges = []

      for sec in range(num_co):
        tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
      print tranges
      time_c = np.mean(tranges, axis=1)
      print 'center time:'
      print time_c

      data = get_data(tranges, scan_name)
      data_star = get_data_star(tranges, scan_name)

      for sec in range(num_co):
        if (len(data[sec]) == 0) or (len(data_star[sec]) == 0):
          continue
        data_sec =  np.array(data[sec], dtype='float64')
        #bin_masks = [data_sec[:,3]<16, data_sec[:,3]>=16, data_sec[:,5]<8.5, data_sec[:,5]>=8.5]
        limit_list_xa = [[0,10], [10,20], [20,50]]
        limit_list_q = [[0,6], [6,11], [11,50]]
        bin_masks = []
        for f in [3,5]:
          if f==3:
            limit = limit_list_xa
          elif f==5:
            limit = limit_list_q
          for bin in range(nbin):
            low = limit[bin][0] #np.percentile(data_sec[:,f], float(bin)/float(nbin)*100.)
            up = limit[bin][1] #np.percentile(data_sec[:,f], float(bin+1)/float(nbin)*100.)
            #print (f, low, up)
            bin_masks.append((data_sec[:,f]>=low)&(data_sec[:,f]<=up))
        for bin_num in range(nbin*2):
          bin_mask = bin_masks[bin_num]
          data_tmp = data_sec[bin_mask]
          coo1 = np.array(data_tmp, dtype='float64')[:,-5:-3]
          coo2 = np.array(data_star[sec], dtype='float64')[:,1:3]

          coo1 = get_detector_pos(coo1, detector_size[0])
          coo2 = get_detector_pos(coo2, detector_size[0])
          print coo1.shape
          print coo2.shape

          for i in range(detector_size[0]):
            for j in range(detector_size[1]):
              mask = (coo1[:,0]>=i) & (coo1[:,0]<i+1) & (coo1[:,1]>=j) & (coo1[:,1]<j+1)
              tmp_coo1 = coo1[mask]

              mask = (coo2[:,0]>=i) & (coo2[:,0]<i+1) & (coo2[:,1]>=j) & (coo2[:,1]<j+1)
              tmp_coo2 = coo2[mask]

              bins[bin_num][i*detector_size[1]+j] += get_corr_map(tmp_coo2, tmp_coo1, hist_size)

    for bin_num in range(nbin*2):
      np.save('../data/%s-cal-sec-star/hists-%d-%d.npy'%(scan_name, bin_num, pid), bins[bin_num])


if __name__ == '__main__':

#hpc single, identify outliers, 0.5s, secondary co
  if True:
    name = sys.argv[1]

    output = '../data/{0}-cal-sec-star/hists.npy'.format(name)
    if os.path.isfile(output):
      print('existed')
      #exit()

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    length = co_data.shape[0]
    
    p_num = int(sys.argv[2])
    chunk_len = int(length/p_num)

    try:
      asp_cal = np.load('../data/photon_list/%s-cal-sec_asp_new.npy'%name)
    except IOError:
      T = co_data['T']
      ra = co_data['ra']
      dec = co_data['dec']
      roll = co_data['roll']
      t_new = np.arange((T.shape[0]-1)*200)*0.005+T[0]
      ra_new = np.interp(t_new, T, ra)
      dec_new = np.interp(t_new, T, dec)
      roll_new = np.interp(t_new, T, roll)
      asp_cal = np.array([t_new, ra_new, dec_new, roll_new]).T
      np.save('../data/photon_list/%s-cal-sec_asp_new.npy'%name, asp_cal)

    hist_size = [100,100]
    detector_size = [50,50]
    nbin = 3

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_r_sec, args=(pid, name, 1, asp_cal, start, end, hist_size, detector_size, nbin, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    #run_one_r_sec(pid, name, 1., start, end, return_dict)

    p = Process(target=run_one_r_sec, args=(pid, name, 1, asp_cal, start, end, hist_size, detector_size, nbin, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'

    bin_list = ['xa_0', 'xa_1', 'xa_2', 'q_0', 'q_1', 'q_2']
    #bin_list = []
    for bin_num in range(nbin*2):
      bin_name = bin_list[bin_num]
      hist_list = []
      for i in range(detector_size[0]):
        for j in range(detector_size[1]):
          hist_list.append(np.zeros(hist_size))

      for pid in range(0, p_num):
        hists = np.load('../data/%s-cal-sec-star/hists-%d-%d.npy'%(name, bin_num, pid))
        for i in range(detector_size[0]):
          for j in range(detector_size[1]):
            hist_list[i*detector_size[1]+j] += hists[i*detector_size[1]+j]
      np.save('../data/%s-cal-sec-star/hists-%s.npy'%(name, bin_name), hist_list)

      for pid in range(0, p_num):
        subprocess.call(['rm', '../data/%s-cal-sec-star/hists-%d-%d.npy'%(name, bin_num, pid)])
    subprocess.call(['rm', '-r', '../data/%s-cal-sec-star/split'%(name)])

    '''
    centroids = []
    for i in range(detector_size[0]):
      for j in range(detector_size[1]):
        H = hist_list[i*detector_size[1]+j]    
        t = 1
        print 'threshold:{0}'.format(t) 
        #plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'))
        #plt.show()
        data = H.byteswap(True).newbyteorder()
        data = data.copy(order='C')
        data = data.byteswap(True).newbyteorder()
        c_data = c3.find_centroid(data, t)
        if c_data is not None:
          cy, cx, max_value, flux = c_data
          centroids.append([-0.2+cx*0.4/hist_size[0], -0.2+cy*0.4/hist_size[0]])
        else:
          centroids.append([0, 0])
    np.save('../data/%s-cal-sec-star/centroids.npy'%(name), centroids)
    '''
