import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import csv
import imagetools
import numpy as np
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

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def get_corr_map_(coo1, coo2, skypos, skyrange, sec):
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
    for i in range(len2):
      print(i)
      tmp_co = np.roll(coo2, i, axis=0)[0:len1,:]-coo1
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=0.01,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=0.01,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
  else:
    for i in range(len1):
      print(i)
      tmp_co = coo2-np.roll(coo1, i, axis=0)[0:len2,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=0.01,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=0.01,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)

  return co_rel[1:]

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
    for i in range(len2):
      #print(i)
      tmp_co = np.roll(coo2, i, axis=0)[0:len1,:]-coo1
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
      '''
      if (i+1)%200 == 0:
        writer.writerows(co_rel[1:])
        foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
        H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                   bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
        count += H
        co_rel = np.array([[0, 0]])
      '''
  else:
    for i in range(len1):
      #print(i)
      tmp_co = coo2-np.roll(coo1, i, axis=0)[0:len2,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
      '''
      if (i+1)%200 == 0:
        writer.writerows(co_rel[1:])
        foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
        H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                   bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
        count += H
        co_rel = np.array([[0, 0]])
      '''
  '''
  if co_rel.shape[0]>1:
    #writer.writerows(co_rel[1:])
    foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
    H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                               bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
    count += H
  #csvfile.close()
  '''
  print co_rel.shape
  if co_rel.shape[0]>50:
    centroid = ck.find_centroid(co_rel[1:], bandwidth, 11, bound)
  else:
    return count, np.array([0.0, 0.0])

  return count, centroid

def get_corr_map_fast(coo1, coo2, skypos, skyrange, sec, bound, bandwidth):
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
  print co_rel.shape
  if co_rel.shape[0]>50:
    centroid = ck.find_centroid(co_rel[1:], bandwidth, 11, bound)
  else:
    return count, np.array([0.0, 0.0])

  return count, centroid

#identify outliers
def get_corr_map_o(coo1, coo2, skypos, skyrange, sec, bound, bandwidth):
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
    for i in range(len2):
      #print(i)
      tmp_co = np.roll(coo2, i, axis=0)[0:len1,:]-coo1
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
      '''
      if (i+1)%200 == 0:
        writer.writerows(co_rel[1:])
        foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
        H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                   bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
        count += H
        co_rel = np.array([[0, 0]])
      '''
  else:
    for i in range(len1):
      #print(i)
      tmp_co = coo2-np.roll(coo1, i, axis=0)[0:len2,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
      '''
      if (i+1)%200 == 0:
        writer.writerows(co_rel[1:])
        foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
        H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                   bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
        count += H
        co_rel = np.array([[0, 0]])
      '''
  '''
  if co_rel.shape[0]>1:
    #writer.writerows(co_rel[1:])
    foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
    H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                               bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
    count += H
  #csvfile.close()
  '''
  print co_rel.shape
  if co_rel.shape[0]>50:
    centroid = ck.find_centroid(co_rel[1:], bandwidth, 11, bound)
  else:
    return count, np.array([0.0, 0.0]), 1

  return count, centroid, 0



def get_corr_map_photon(coo1, coo2, skypos, skyrange, sec, bound, bandwidth):
  imsz = imagetools.deg2pix(skypos, skyrange, 0.0001)
  count = np.zeros(imsz)
  print(imsz)
  co_rel = np.array([[0,0]])
  len1 = coo1.shape[0]
  len2 = coo2.shape[0]
  print(len1,len2)
  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

  if len2<=80 or len1<=80:
    return count, np.array([0.0, 0.0])
  if len2>len1:
    for i in range(len2):
      #print(i)
      tmp_co = np.roll(coo2, i, axis=0)[0:len1,:]-coo1
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)

  else:
    for i in range(len1):
      #print(i)
      tmp_co = coo2-np.roll(coo1, i, axis=0)[0:len2,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)

  print co_rel.shape
  if co_rel.shape[0]>1000:
    centroid = ck.find_centroid(co_rel[1:], bandwidth, 11, bound)
  else:
    return count, np.array([0.0, 0.0])

  return count, centroid

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

def run_one(pid, scan_name, step, return_dict):

    num_co = int(1/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    print length
    for initial_sec in range(1, length):
      file_path = '../data/%s/cata/centroids%d.npy'%(scan_name, initial_sec)
      if os.path.exists(file_path):
        print 'skip:%d'%initial_sec
        continue

      intitial_asp = co_data[initial_sec]
      center = np.array([intitial_asp[1], intitial_asp[2]])  

      print(intitial_asp)

      skypos = [0.0, 0.0]
      skyrange = [0.2, 0.2]

      wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

      initial_time = intitial_asp[0]-0.5

      tranges = []

      for sec in range(num_co):
        tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
      print tranges

      data = get_data(tranges, scan_name)
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

        count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.1, 0.008)  
        centroids.append(centroid)
        print centroid
      print centroids
      np.save('../data/%s/cata/centroids%d.npy'%(scan_name, initial_sec), centroids)

def run_one_r(pid, scan_name, step, start, end, return_dict):

    num_co = int(1/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    print length
    for initial_sec in range(start, end):
      file_path = '../data/%s/cata/centroids%d_half.npy'%(scan_name, initial_sec)
      if os.path.exists(file_path):
        print 'skip:%d'%initial_sec
        continue

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

      data = get_data(tranges, scan_name)
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

        count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.15, 0.008)  
        centroids.append(centroid)
        print centroid
      print centroids
      np.save('../data/%s/cata/centroids%d_half.npy'%(scan_name, initial_sec), centroids)

def run_one_r_sec(pid, scan_name, step, start, end, return_dict):

    num_co = int(1/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    cal_initial = int((co_data[1][0]-0.5)*1000)
    offsets = np.load('../data/%s/cata/offsets1_10_new_inter_half_fine.npy'%scan_name)

    print length
    for initial_sec in range(start, end):
      file_path = '../data/%s/cata/centroids_sec%d.npy'%(scan_name, initial_sec)
      if os.path.exists(file_path):
        print 'skip:%d'%initial_sec
        continue

      intitial_asp = co_data[initial_sec]
      #center = np.array([intitial_asp[1], intitial_asp[2]])
      #use the first order calibrated asp instead of the original one
      center = cal_asp_r(offsets, cal_initial, intitial_asp[0]*1000, intitial_asp[1:3], 200) 

      print(intitial_asp)

      skypos = [0.0, 0.0]
      skyrange = [0.1, 0.1]

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

        count, centroid = get_corr_map_fast(coo2, coo1, skypos, skyrange, sec, 0.03, 0.0015)  
        centroids.append(centroid)
        print centroid
      print centroids
      np.save('../data/%s/cata/centroids_sec%d.npy'%(scan_name, initial_sec), centroids)


#increase grid size, identify outliers
def run_one_r_o(pid, scan_name, step, start, end, return_dict):

    num_co = int(1/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    print length
    for initial_sec in range(start, end):
      file_path = '../data/%s/cata/centroids_new%d.npy'%(scan_name, initial_sec)
      if os.path.exists(file_path):
        print 'skip:%d'%initial_sec
        continue

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

      data = get_data(tranges, scan_name)
      print len(data)

      output = "../data/%s/cata/%s.csv"%(scan_name, scan_name)
      dir = os.path.dirname(output)
      if not os.path.exists(dir):
          os.makedirs(dir)
      centroids = []
      out_mask = []
      for sec in range(num_co):
        if len(data[sec]) ==0:
          centroid = np.array([0.0, 0.0])
          out_mask.append(1)
          centroids.append(centroid)
          print centroid
          continue
        coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]
        coo2 = catalog_fits.get_catalog(center, 0.69)

        coo1 = angle_filter(coo1, center, 1.)
        coo2 = angle_filter(coo2, center, 1.)

        count, centroid, out = get_corr_map_o(coo2, coo1, skypos, skyrange, sec, 0.15, 0.008)  
        centroids.append(centroid)
        out_mask.append(out)
        print centroid
      print centroids
      np.save('../data/%s/cata/centroids_new%d.npy'%(scan_name, initial_sec), centroids)
      np.save('../data/%s/cata/out_mask%d.npy'%(scan_name, initial_sec), out_mask)

def get_data_star(scan_name, star):
  data = []
  time_list = []
  count = []
  csv_file = '../data/%s/star/list/%s.csv'%(scan_name, star)
  with open(csv_file, 'rb') as file:
    reader = csv.reader(file)
    first = reader.next()
    last = int(float(first[0])/1000.)
    time_list.append(last)
    data.append([])
    data[0].append(first)
    i = 0
    count.append(1)
    for row in reader:
      time = int(float(row[0])/1000.)
      if time - last > 0:
        data.append([])
        i += 1
        data[i].append(row)
        last = time
        time_list.append(last)
        count.append(1)
      else:
        data[i].append(row)
        count[-1] += 1

  return data, time_list ,count

def run_one_photon(pid, scan_name, start, end, return_dict):

  cata = spi.load_obj('../data/%s_starset_new'%scan_name)
  cata_a = np.array(list(cata))
  cata_len = len(cata)

  skypos = [0.0, 0.0]
  skyrange = [0.2, 0.2]

  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

  for star_num in range(start, end+1): 
    order = 1
    star = 'star_%d-%d'%(star_num, order)

    if not os.path.exists('../data/%s/star/list/%s.csv'%(scan_name, star)):
      print 'skip:%s'%star
      continue

    data, time_list, count = get_data_star(scan_name, star)
    length =  len(data)
    mean = np.mean(count)
    if mean <= 80:
      print 'skip:%s'%star
      continue

    output = "../data/%s/star/co_photon/%s.csv"%(scan_name, scan_name)
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
        os.makedirs(dir)
    centroids = []
    for sec in range(length):
      if len(data[sec]) ==0:
        centroid = np.array([0.0, 0.0])
        centroids.append(centroid)
        print centroid
        continue
      coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]
      #coo2 = np.array(data[sec+1], dtype='float64')[:,-3:-1]
      coo2 = np.array([cata_a[star_num,:]])

      count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.03, 0.004)  
      centroids.append(centroid)
      print centroid
    print centroids
    np.save('../data/%s/star/co_photon/centroids%d-%d_c.npy'%(scan_name, star_num, order), centroids)

    time = np.array(time_list)
    centers = np.array(centroids)

    '''
    for i in range(1,centers.shape[0]):
      centers[i,:] += centers[i-1,:]
    '''

    print time
    print centers
    f, axes = plt.subplots(2, 1)
    axes[0].plot(time, centers[:,0])
    axes[0].set_ylabel(r'$\delta RA/degree$')
    ylim = axes[0].get_ylim()
    axes[0].set_ylim(ylim+np.absolute(ylim)*0.01)
    plt.setp( axes[0].get_xticklabels(), visible=False)
    #plt.setp( axes[0].get_yticklabels(), visible=False)
    
    axes[1].plot(time, centers[:,1])
    axes[1].set_ylabel(r'$\delta DEC/degree$')
    axes[1].set_xlabel('time/s')

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0, hspace=0)
    plt.savefig('../plots/%s/star/co_photon/offset%d-%d_c.png'%(scan_name,star_num,order),dpi=190)  
    plt.clf()
    '''
    plt.plot(time, centers[:,0])
    plt.xlabel('time/s')
    plt.ylabel(r'$\delta RA/degree$')
    plt.savefig('../plots/%s/star/co_photon/ra%d-%d_c.png'%(scan_name,star_num,order),dpi=190)
    plt.clf()
    
    plt.plot(time, centers[:,1])
    plt.xlabel('time/s')
    plt.ylabel(r'$\delta DEC/degree$')
    plt.savefig('../plots/%s/star/co_photon/dec%d-%d_c.png'%(scan_name,star_num,order),dpi=190)
    plt.clf()
    '''

def run_one_photon_half(pid, scan_name, step, start, end, return_dict):

    num_co = int(0.5/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    print length
    for initial_sec in range(start, end):
      file_path = '../data/%s/cata/centroids_half_photon%d.npy'%(scan_name, initial_sec)
      if os.path.exists(file_path):
        print 'skip:%d'%initial_sec
        continue 

      intitial_asp = co_data[initial_sec]
      print(intitial_asp)

      center = np.array([intitial_asp[1], intitial_asp[2]])
      centroids = []
      for part in range(2):
        print 'part:%d'%part
        skypos = [0.0, 0.0]
        skyrange = [0.1, 0.1]

        wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

        initial_time = intitial_asp[0]-0.5

        tranges = []

        for sec in range(num_co):
          tranges.append([initial_time+step*sec+part*0.5, initial_time+step*(sec+1)+part*0.5])
        print tranges

        data = get_data(tranges, scan_name)
        print len(data)

        output = "../data/%s/cata/%s.csv"%(scan_name, scan_name)
        dir = os.path.dirname(output)
        if not os.path.exists(dir):
            os.makedirs(dir)

        if len(data[2]) == 0:
          for sec in range(4):
            centroid = np.array([0.0, 0.0])
            centroids.append(centroid)
            print centroid
          continue
        else:
          coo2 = np.array(data[2], dtype='float64')[:,-3:-1]
          coo2 = angle_filter(coo2, center, 0.4)

        for sec in [0,1,3,4]:
          if len(data[sec]) ==0:
            centroid = np.array([0.0, 0.0])
            centroids.append(centroid)
            print centroid
            continue
          coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]

          coo1 = angle_filter(coo1, center, 0.4)

          count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.05, 0.004)  
          centroids.append(centroid)
          print centroid
      print centroids
      np.save('../data/%s/cata/centroids_half_photon%d.npy'%(scan_name, initial_sec), centroids)

if __name__ == '__main__':
  if False:
    initial_sec = int(sys.argv[1])

    if len(sys.argv)>2:
      step = float(sys.argv[2])
    else:
      step = 0.5

    num_co = int(1/step)

    hdulist = pyfits.open('../data/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data
    intitial_asp = co_data[initial_sec]
    center = np.array([intitial_asp[1], intitial_asp[2]])  

    print(intitial_asp)

    skypos = [0.0, 0.0]
    skyrange = [0.04, 0.04]

    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

    initial_time = intitial_asp[0]-0.5

    tranges = []

    for sec in range(num_co+1):
      tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
    print tranges

    data = get_data(tranges)
    print len(data)

    centroids = []
    for sec in range(num_co):
      coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]
      coo2 = np.array(data[sec+1], dtype='float64')[:,-3:-1]

      coo1 = angle_filter(coo1, center, 0.4)
      coo2 = angle_filter(coo2, center, 0.4)

      count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.02, 0.001)  
      centroids.append(centroid)
      print centroid

    np.save('../data/05_r/centroids%d.npy'%(initial_sec), centroids) 

  if False:
    initial_sec = int(sys.argv[1])

    if len(sys.argv)>3:
      step = float(sys.argv[2])
    else:
      step = 0.5

    scan_name = sys.argv[-1]

    num_co = int(1/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    intitial_asp = co_data[initial_sec]
    center = np.array([intitial_asp[1], intitial_asp[2]])  

    print(intitial_asp)

    skypos = [0.0, 0.0]
    skyrange = [0.2, 0.2]

    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

    initial_time = intitial_asp[0]-0.5

    tranges = []

    for sec in range(num_co):
      tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
    print tranges

    data = get_data(tranges, scan_name)
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

      count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.1, 0.008)  
      centroids.append(centroid)
      print centroid
    print centroids

    np.save('../data/%s/cata/centroids%d.npy'%(scan_name, initial_sec), centroids)

#for hpc
  if False:
    initial_sec = int(sys.argv[1])

    if len(sys.argv)>3:
      step = float(sys.argv[2])
    else:
      step = 0.5

    scan_name = sys.argv[-1]

    num_co = int(1/step)

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    print length
    for initial_sec in range(1, length):
      file_path = '../data/%s/cata/centroids%d.npy'%(scan_name, initial_sec)
      if os.path.exists(file_path):
        print 'skip:%d'%initial_sec
        continue

      intitial_asp = co_data[initial_sec]
      center = np.array([intitial_asp[1], intitial_asp[2]])  

      print(intitial_asp)

      skypos = [0.0, 0.0]
      skyrange = [0.2, 0.2]

      wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0001)

      initial_time = intitial_asp[0]-0.5

      tranges = []

      for sec in range(num_co):
        tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
      print tranges

      data = get_data(tranges, scan_name)
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

        count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.1, 0.008)  
        centroids.append(centroid)
        print centroid
      print centroids

      np.save('../data/%s/cata/centroids%d.npy'%(scan_name, initial_sec), centroids)

#for hpc process
  if False:
    name_file = sys.argv[1]
    with open('../name/%s'%name_file) as f:
      name_list = f.read().splitlines()
    print name_list

    p_num = len(name_list)

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for name in name_list:
      p = Process(target=run_one, args=(pid, name, 0.1, return_dict))
      p.start()
      p_list.append(p)
      pid += 1

    for p in p_list:
      p.join()
    print 'all done'

#hpc single
  if False:
    name = sys.argv[1]

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    length = co_data.shape[0]
    
    p_num = 10
    chunk_len = int(length/p_num)

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_r, args=(pid, name, 0.1, start, end, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    p = Process(target=run_one_r, args=(pid, name, 0.1, start, end, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'

#hpc single, increase grid size, identify outliers
  if False:
    name = sys.argv[1]

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    length = co_data.shape[0]
    
    p_num = 10
    chunk_len = int(length/p_num)

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_r_o, args=(pid, name, 0.1, start, end, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    p = Process(target=run_one_r_o, args=(pid, name, 0.1, start, end, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'

#hpc single, identify outliers, 0.5s
  if False:
    name = sys.argv[1]

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    length = co_data.shape[0]
    
    p_num = 11
    chunk_len = int(length/p_num)

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_r, args=(pid, name, 0.5, start, end, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    p = Process(target=run_one_r, args=(pid, name, 0.5, start, end, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'

#hpc single, identify outliers, 0.5s, secondary co
  if True:
    name = sys.argv[1]

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    length = co_data.shape[0]
    
    p_num = int(sys.argv[2])
    chunk_len = int(length/p_num)

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_r_sec, args=(pid, name, 1., start, end, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    p = Process(target=run_one_r_sec, args=(pid, name, 1., start, end, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'
    centroid_csv.generate_sec_offsets(name)

