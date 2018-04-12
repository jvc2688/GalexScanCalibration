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
import centroid_kde_d as ck
import pickle

def load_obj(name ):
    with open('../data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def get_corr_map(coo1, coo2, skypos, skyrange, sec):
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
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=0.03,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=0.03,:]
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
      print(i)
      tmp_co = coo2-np.roll(coo1, i, axis=0)[0:len2,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=0.03,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=0.03,:]
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
  centroid = ck.find_centroid(co_rel[1:], 0.001, 12)

  return count, centroid

def get_data(tranges):
  look_up = load_obj('lookup')
  data = []
  trange = [0 , 0]
  num_interval = len(tranges)
  tranges = np.array(tranges, dtype='double')
  tranges = (tranges*1000).astype(int)
  initial_time = trange[0] = int(tranges[0,0]) #int(reader.next()[0]) + 0
  final_time = trange[1] = int(tranges[num_interval-1, 1])
  start_index = look_up['%d'%initial_time]
  end_index = look_up['%d'%final_time]
  csv_list = []
  for i in range(start_index, end_index+1):
    csv_list.append('../data/split/%d.csv'%i)
 
  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)
      for i in range(num_interval):
        data.append([])
      print len(data)

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

    coo1 = angle_filter(coo1, center, 1.)
    coo2 = angle_filter(coo2, center, 1.)

    count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec)  
    centroids.append(centroid)
    print centroid

  np.save('../data/centroids%d.npy'%(initial_sec), centroids) 
