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
import centroid_plot_csv as centroid_csv
import gnomonic as gn
import c3
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude


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

def get_corr_map(coo1, coo2, skypos, skyrange, sec, bound, bandwidth):
  imsz = imagetools.deg2pix(skypos, skyrange, 0.0001)
  count = np.zeros(imsz)
  print(imsz)
  co_rel = np.array([[0,0]])
  len1 = coo1.shape[0]
  len2 = coo2.shape[0]
  print(len1,len2)
  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)
  #with open('../data/try2_%d.csv'%sec, 'wb') as csvfile:
    #writer = csv.writer(csvfile)
  if len2>len1:
    for i in range(len1):
      #print(i)
      tmp_co = coo2-coo1[i,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
  else:
    for i in range(len2):
      #print(i)
      tmp_co = coo2[i,:]-coo1
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)

  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                             bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))


  data = H.byteswap(True).newbyteorder()
  data = data.copy(order='C')
  data = data.byteswap(True).newbyteorder()
  cy, cx, max_value = c3.find_centroid(data)
  centroid = wcs.wcs_pix2world(wcs.sip_foc2pix([[cx, cy]],1),1)[0]
  if centroid[0]>1:
    centroid[0] = centroid[0]-360.

  return centroid, max_value


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
  look_up = load_obj('../data/%s/split/lookup'%scan_name)
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
    csv_list.append('../data/%s/split/%d.csv'%(scan_name, i))
 
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

def get_data_cal(tranges, scan_name, cal_initial, offsets):

  step = 200

  look_up = load_obj('../data/%s/split/lookup'%scan_name)
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
    csv_list.append('../data/%s/split/%d.csv'%(scan_name, i))
 
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
            row = cal_photon_r(offsets, cal_initial, time, row, step)
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

def run_one_r(pid, scan_name, step, start, end, return_dict):

    print 'run_one_r'
    num_co = int(1/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    plt.plot(co_data['T'])
    plt.show()

    cal_initial = int((co_data[1][0]-0.5)*1000)
    print cal_initial
    offsets = np.load('../data/%s/cata/offsets1_10_new_inter_half_fine.npy'%scan_name)

    print length
    print co_data[-10]
    print co_data[0]
    print co_data[-1][0]-co_data[0][0]
    #print offsets[300*200:301*200]
    for initial_sec in range(start, end):
      file_path = '../data/%s/cata/centroids%d_half.npy'%(scan_name, initial_sec)
      '''
      if os.path.exists(file_path):
        print 'skip:%d'%initial_sec
        continue
      '''
      intitial_asp = co_data[initial_sec]
      center = np.array([intitial_asp[1], intitial_asp[2]])  

      print(intitial_asp)

      skypos = [0.0, 0.0]
      skyrange = [0.3, 0.3]

      wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

      initial_time = intitial_asp[0]-0.5

      tranges = []

      for sec in range(num_co):
        tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
      print tranges

      data = get_data_cal(tranges, scan_name, cal_initial, offsets)
      #data = get_data(tranges, scan_name)
      print len(data)

      output = "../data/%s/cata/%s.csv"%(scan_name, scan_name)
      dir = os.path.dirname(output)
      if not os.path.exists(dir):
          os.makedirs(dir)
      centroids = []
      for sec in range(num_co):
        if len(data[sec]) ==0:
          centroid = np.array([0.0, 0.0])
          centroids.append(centroid)
          print centroid
          continue
        coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]
        coo2 = catalog_fits.get_catalog(center, 0.69)

        coo1 = angle_filter(coo1, center, 1.)
        coo2 = angle_filter(coo2, center, 1.)

        get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.05, 0.008)
      '''
        count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.15, 0.008)  
        centroids.append(centroid)
        print centroid
      print centroids
      np.save('../data/%s/cata/centroids%d_half.npy'%(scan_name, initial_sec), centroids)
      '''

def run_one_r_sec(pid, scan_name, step, start, end, return_dict):

    pixsz = 0.000416666666666667
    skypos_count = [5.,-9.5]
    skyrange_count = [2., 3.]
    tranges_count = [[0,0]]
    imsz_count = imsz = imagetools.deg2pix_g(skypos_count, skyrange_count, pixsz)
    wcs_count = imagetools.define_wcs_g(skypos_count,skyrange_count,width=False,height=False,verbose=0,pixsz=pixsz)
    count = np.zeros(imsz_count)

    num_co = int(1/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    cal_initial = int((co_data[1][0]-0.5)*1000)
    offsets = np.load('../data/%s/cata/offsets_inter_half_sec.npy'%scan_name)

    angle_list = [-0.035, -0.03, -0.025, -0.02, -0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035]

    print length
    for initial_sec in range(start, end):
      print 'time:%d'%initial_sec
      file_path = '../data/%s/cata/centroids_sec%d.npy'%(scan_name, initial_sec)
      '''
      if os.path.exists(file_path):
        print 'skip:%d'%initial_sec
        continue
      '''
      intitial_asp = co_data[initial_sec]
      rot = intitial_asp[3]
      center = np.array([intitial_asp[1], intitial_asp[2]])
      #use the first order calibrated asp instead of the original one
      center = cal_asp_r(offsets, cal_initial, intitial_asp[0]*1000, intitial_asp[1:3], 200) 

      #print(intitial_asp)

      skypos = [0.0, 0.0]
      skyrange = [0.02, 0.02]

      wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

      initial_time = intitial_asp[0]-0.5

      tranges = []

      for sec in range(num_co):
        tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
      print tranges

      data = get_data_cal(tranges, scan_name, cal_initial, offsets)
      print len(data)

      output = "../data/%s/cata/%s.csv"%(scan_name, scan_name)
      dir = os.path.dirname(output)
      if not os.path.exists(dir):
          os.makedirs(dir)
      centroids = []
      for sec in range(num_co):
        if len(data[sec]) ==0:
          centroid = np.array([0.0, 0.0])
          centroids.append(centroid)
          print centroid
          continue
        coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]
        coo2 = catalog_fits.get_catalog(center, 0.69)

        coo1 = angle_filter(coo1, center, 1.)
        coo2 = angle_filter(coo2, center, 1.)


        xi, eta = gn.gnomfwd_simple(coo1[:,0], coo1[:,1], center[0], center[1], -rot, 1/36000., 0.0)

        max_now = -1000
        cent_now = []
        for angle in angle_list:
          coo1[:,0], coo1[:,1] = gn.gnomrev_simple(xi,eta,center[0],center[1],-(rot+angle),1/36000.,0.)
          centroid, max_value = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.02, 0.004)
          if max_value>max_now:
            max_now = max_value
            cent_now = np.append(centroid, angle) #centroid.append(angle)
        

        coo1[:,0], coo1[:,1] = gn.gnomrev_simple(xi,eta,center[0]-cent_now[0],center[1]-cent_now[1],-(rot+cent_now[2]),1/36000.,0.)
        
        sky_data = SkyCoord(coo1, unit='deg', frame=FK5, equinox='J2000.0')
        gal = sky_data.transform_to(Galactic)
        data_count = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1) 
        H = hist_g(data_count, wcs_count, imsz_count)
        count += H
        sky_data=gal=data_count=H = None
        gc.collect()

        #xi2, eta2 = gn.gnomfwd_simple(coo2[:,0], coo2[:,1], center[0], center[1], -rot, 1/36000., 0.0)

        #r, phi = cart2pol(xi, eta)
        #r2, phi2 = cart2pol(xi2, eta2)

        #coo1 = np.array([r,phi]).T
        #coo2 = np.array([r2,phi2]).T

        #get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.05, 0.0015)  

        centroids.append(cent_now)
        print cent_now
      print centroids
      #np.save('../data/%s/cata/centroids_rot%d.npy'%(scan_name, initial_sec), centroids)
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos_count, tranges_count, skyrange_count, hdu=hdu, wcs=wcs_count)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_%s_gal_sec_rot_count.fits'%(scan_name), clobber=False)


if __name__ == '__main__':

#hpc single, identify outliers, 0.5s, secondary co
  if True:
    name = sys.argv[1]

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    length = co_data.shape[0]
    
    p_num = 1
    chunk_len = int(length/p_num)

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_sec, args=(pid, name, 1., start, end, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    run_one_r_sec(pid, name, 1., start, end, return_dict)
    '''
    p = Process(target=run_one_r_sec, args=(pid, name, 1., start, end, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'
    centroid_csv.generate_sec_offsets(name)
    '''

