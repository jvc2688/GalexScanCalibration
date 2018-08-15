import matplotlib
matplotlib.use('Agg')
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
from scipy import interpolate


def dis_correct(data, dis_map, asp):
    ix = np.digitize(data[:,0]/1000., asp[:,0])-1

    low_mask = data[:,3]<16
    high_mask = data[:,3]>=16

    coo = data[low_mask, 6:8]
    data[low_mask, 6] -= dis_map['xi_l'](coo)
    data[low_mask, 7] -= dis_map['eta_l'](coo)

    coo = data[high_mask, 6:8]
    data[high_mask, 6] -= dis_map['xi_h'](coo)
    data[high_mask, 7] -= dis_map['eta_h'](coo)

    ra, dec = gn.gnomrev_simple(data[:,6], data[:,7], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return np.array([ra,dec]).T

def dis_correct_non(data, dis_map, asp):
    ix = np.digitize(data[:,0]/1000., asp[:,0])-1

    ra, dec = gn.gnomrev_simple(data[:,6], data[:,7], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return np.array([ra,dec]).T

def no_dis(data, dis_map, asp):
    ix = np.digitize(data[:,0]/1000., asp[:,0])-1

    low_mask = data[:,3]<16
    high_mask = data[:,3]>=16

    coo = data[low_mask, 6:8]
    data[low_mask, 6] -= dis_map['xi_l'](coo)
    data[low_mask, 7] -= dis_map['eta_l'](coo)

    coo = data[high_mask, 6:8]
    data[high_mask, 6] -= dis_map['xi_h'](coo)
    data[high_mask, 7] -= dis_map['eta_h'](coo)

    ra, dec = gn.gnomrev_simple(data[:,6], data[:,7], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return np.array([ra,dec]).T

def detector2gon(pos_list, size):
  pos_array = np.array(pos_list, dtype='float64')
  pos = (pos_array/50.-0.5)*36000*800*0.001666
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

def get_corr_map(coo1, coo2, skypos, skyrange, sec, pixsz, idx, suffix):
  imsz = imagetools.deg2pix(skypos, skyrange, pixsz).astype(int)
  count = np.zeros(imsz)
  #print(imsz)
  bound = skyrange[0]
  co_rel = np.array([[0,0]])
  len1 = coo1.shape[0]
  len2 = coo2.shape[0]
  print(len1,len2)
  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
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

  t = np.percentile(H, 90)
  t = 1
  print 'threshold:{0}'.format(t) 
  hp = H.copy()

  #plt.colorbar()
  
  #plt.plot(co_rel[:,0],co_rel[:,1], 'o', alpha=0.05)
  #plt.xlim(-0.005,0.005)
  #plt.ylim(-0.005,0.005)

  data = H.byteswap(True).newbyteorder()
  data = data.copy(order='C')
  data = data.byteswap(True).newbyteorder()
  c_data = c3.find_centroid(data, t)
  if c_data is not None:
    cy, cx, max_value, flux = c_data
    cy+=0.5
    cx+=0.5
    centroid = wcs.wcs_pix2world(wcs.sip_foc2pix([[cx, cy]],1),1)[0]
  else:
    centroid = [0.,0.]
    max_value = 0
    flux = 500
  if centroid[0]>1:
    centroid[0] = centroid[0]-360.
  #print 'max:{0}'.format(max_value)
  print hp.shape

  plt.imshow(hp,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[0.02, -0.02, 0.02, -0.02], origin='upper', imlim=[-0.005,0.005, 0.02, -0.02])
  plt.colorbar()
  plt.xlim(-0.005, 0.005)
  plt.ylim(-0.005, 0.005)
  plt.ylabel('gb')
  plt.xlabel('gl')
  plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
  plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
  plt.plot(centroid[0],centroid[1],'+r', markersize=8)

  #plt.show()
  plt.savefig('/home/dw1519/galex/plots/co14non_c/{0}.pdf'.format(idx), dpi=190)
  plt.clf()

  return centroid, max_value, flux


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
      skyrange = [0.2, 0.2]

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

        coo1 = angle_filter(coo1, center, 0.5)
        coo2 = angle_filter(coo2, center, 0.5)

        get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.05, 0.008)
      '''
        count, centroid = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.15, 0.008)  
        centroids.append(centroid)
        print centroid
      print centroids
      np.save('../data/%s/cata/centroids%d_half.npy'%(scan_name, initial_sec), centroids)
      '''

def run_one_r_sec(pid, scan_name, step, resolution, asp_cal, start, end, dis_map, return_dict):

    print('run one r sec')

    num_co = int(resolution/step)
    
    #load asp file
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    co_data = hdulist[1].data
    length = co_data.shape[0]

    '''
    try:
      asp_cal = np.load('../data/photon_list/%s_asp_new.npy'%scan_name)
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
      np.save('../data/photon_list/%s_asp_new.npy'%scan_name, asp_cal)
    '''

    #cal_initial = int((co_data[1][0]-0.5)*1000)
    #offsets = np.load('../data/%s/cata/offsets_inter_half_sec.npy'%scan_name)

    angle_list = [0.]

    print length
    data1 = []
    data2 = []
    start = 700
    c_list = []
    move = 100
    for start in range(1, 1200, move):
      print start
      data1 = []
      data2 = []
      for initial_sec in range(start, start+move):
        print 'time:%d'%initial_sec
        file_path = '../data/%s/cata/centroids_rot%d.npy'%(scan_name, initial_sec)

        intitial_asp = co_data[initial_sec]
        rot = intitial_asp[3]
        #center = np.array([intitial_asp[1], intitial_asp[2]])
        #use the first order calibrated asp instead of the original one
        #center = cal_asp_r(offsets, cal_initial, intitial_asp[0]*1000, intitial_asp[1:3], 200) 
        #print(intitial_asp)

        skypos = [0.0, 0.0]
        skyrange = [0.02, 0.02]
        if step == 1:
          skyrange = [0.04, 0.04]#[0.3, 0.3]
        elif step == 0.5:
          skyrange = [0.02, 0.02]

        initial_time = intitial_asp[0]-step/2.

        tranges = []

        for sec in range(num_co):
          tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
        print tranges
        time_c = np.mean(tranges, axis=1)
        print 'center time:'
        print time_c

        #data = get_data_cal(tranges, scan_name, cal_initial, offsets)
        data = get_data(tranges, scan_name)
        print 'number of the photons:'
        print len(data)

        centroids = []
        center_time = []
        for sec in range(num_co):
          center_time.append(time_c[sec])
          if len(data[sec]) == 0:
            centroid = np.array([0.0, 0.0, 0.0])
            centroids.append(centroid)
            print centroid
            continue
          arg = np.argmin(np.absolute(asp_cal[:,0]-time_c[sec]))
          center = asp_cal[arg, 1:3]

          data_sec = dis_correct_non(np.array(data[sec], dtype='float64'), dis_map, asp_cal)
          coo1 = np.array(data_sec, dtype='float64')
          aperture = 0.69
          #coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]
          coo2 = catalog_fits.get_catalog_tycho(center, aperture)

          coo1 = angle_filter(coo1, center, 0.6)
          coo2 = angle_filter(coo2, center, 0.6)

          '''
          fig = plt.gcf()
          fig.set_size_inches(8,8)
          plt.plot(coo1[:,0], coo1[:,1], '.k', markersize=0.1)
          plt.plot(coo3[:,0], coo3[:,1], '+r', markersize=12)
          plt.plot(coo2[:,0], coo2[:,1], '+b', markersize=8)
          #plt.show()
          plt.ylabel('Dec')
          plt.xlabel('RA')
          plt.xlim(center[0]-0.65, center[0]+0.65)
          plt.ylim(center[1]-0.65, center[1]+0.65)
          plt.tight_layout()
          plt.savefig('/home/dw1519/galex/plots/photons.png', dpi=190)
          plt.clf()
          plt.close()
          '''
          '''
          #plot field
          #plt.plot(asp_cal[:,1],'.b')
          #plt.show()
          print center
          print coo1.shape
          print np.min(coo1[:,0]), np.max(coo1[:,0])
          print np.min(coo1[:,1]), np.max(coo1[:,1])
          wcs = imagetools.define_wcs(center,[1.5,1.5],width=False,height=False,verbose=0,pixsz=0.002)
          imsz = imagetools.deg2pix(center, [1.5,1.5], 0.002)
          foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo1,1),1)
          H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                     bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))

          H_new = H.copy()
          data = H_new.byteswap(True).newbyteorder()
          data = data.copy(order='C')
          data = data.byteswap(True).newbyteorder()
          c_data = c3.find_source(data, 1)
          if c_data is not None:
            cx, cy, max_value = c_data
            print cy, cx
            centroid = wcs.wcs_pix2world(wcs.sip_foc2pix([[cx, cy]],1),1)[0]


          coo1 = angle_filter_out(coo1, centroid, 0.01)
          foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo1,1),1)
          H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                     bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))

          plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), vmax=10)
          plt.colorbar()
          plt.show()


          foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo2,1),1)
          H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                     bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
          plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'))
          plt.show()
          '''

          #xi, eta = gn.gnomfwd_simple(coo1[:,0], coo1[:,1], center[0], center[1], -rot, 1/36000., 0.0)

          '''
          max_now = -1000
          cent_now = []
          for angle in angle_list:
            centroid, max_value, flux = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.0002, initial_sec, '1')
            if max_value>max_now:
              max_now = max_value
              cent_now = np.append(centroid, angle) #centroid.append(angle)

          centroids.append(cent_now)
        print centroids
        '''
        #np.save('../data/%s/cata/time_rot%d.npy'%(scan_name, initial_sec), center_time)
        #np.save('../data/%s/cata/centroids_rot%d.npy'%(scan_name, initial_sec), centroids)
        data1.append(coo1)
        data2.append(coo2)
      data1 = np.concatenate(data1, axis=0)
      data2 = np.concatenate(data2, axis=0)
      ar, index = np.unique(data2[:,0], return_index=True)
      data2 = data2[index]

      #sky_data = SkyCoord(data1, unit='deg', frame=FK5, equinox='J2000.0')
      #gal = sky_data.transform_to(Galactic)
      #data1 = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)

      #sky_data = SkyCoord(data2, unit='deg', frame=FK5, equinox='J2000.0')
      #gal = sky_data.transform_to(Galactic)
      #data2 = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)


      centroid, max_value, flux = get_corr_map(data2, data1, skypos, skyrange, sec, 0.0002, initial_sec, '{0}'.format(start))
      print centroid
      c_list.append(centroid)
    np.save('/home/dw1519/galex/plots/co14non_c/c_list.npy', c_list)

if __name__ == '__main__':

#hpc single, identify outliers, 0.5s, secondary co
  if True:
    name = sys.argv[1]

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    co_data = hdulist[1].data
    length = co_data.shape[0]
    
    p_num = int(sys.argv[2])
    chunk_len = int(length/p_num)

    step = float(sys.argv[3])
    resolution = float(sys.argv[4])

    try:
      asp_cal = np.load('../data/photon_list/%s-new_asp.npy'%name)
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
      np.save('../data/photon_list/%s-new_asp_new.npy'%name, asp_cal)

    output = "../data/%s/cata/%s.csv"%(name, name)
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
        os.makedirs(dir)

    detector_size = [50,50]

    dis_map_l = np.load('../data/distortion/194-437-xa_low-centroids.npy')
    dis_map_h = np.load('../data/distortion/194-437-xa_high-centroids.npy')

    #dis_map_l = np.load('../data/distortion/5-185-xa_low-centroids.npy')
    #dis_map_h = np.load('../data/distortion/5-185-xa_high-centroids.npy')

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


    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_r_sec, args=(pid, name, step, resolution, asp_cal, start, end, dis_map, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    #run_one_r_sec(pid, name, 1., start, end, return_dict)

    p = Process(target=run_one_r_sec, args=(pid, name, step, resolution, asp_cal, start, end, dis_map, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'

    #centroid_csv.generate_first_offsets(name)
    #centroid_csv.generate_new_offsets_new(name)

