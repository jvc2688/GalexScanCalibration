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
from scipy import interpolate
from matplotlib.colors import LogNorm
from scipy.interpolate import splev, splrep
from matplotlib.patches import Ellipse
import argparse
import re



def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def dis_correct(data, asp, dis_map=None, slope=None, inter=None, cut=False):
    if cut:
      ya_mask = data[:,4]>=2
      q_mask = data[:,5]>5
      data = data[ya_mask&q_mask]
    ix = np.digitize(data[:,0]/1000., asp[:,0])-1

    if dis_map is not None:
      low_mask = data[:,3]<16
      high_mask = data[:,3]>=16

      coo = data[low_mask, 6:8]
      data[low_mask, 6] -= dis_map['xi_l'](coo)
      data[low_mask, 7] -= dis_map['eta_l'](coo)

      coo = data[high_mask, 6:8]
      data[high_mask, 6] -= dis_map['xi_h'](coo)
      data[high_mask, 7] -= dis_map['eta_h'](coo)

    if slope is not None:
      data[:,7] -= (slope*data[:,4]+inter)*36000*800*0.001666/2400.

    ra, dec = gn.gnomrev_simple(data[:,6], data[:,7], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return np.array([ra,dec]).T, data[:,0]

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

def get_corr_map(coo1, coo2, skypos, skyrange, sec, pixsz, idx, suffix, weight, flux=None):
  imsz = imagetools.deg2pix(skypos, skyrange, pixsz).astype(int)
  count = np.zeros(imsz)
  #print(imsz)
  weights = np.array([])
  bound = skyrange[0]/2.
  co_rel = np.array([[0,0]])
  len1 = coo1.shape[0]
  len2 = coo2.shape[0]
  print(len1,len2)
  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
  #with open('../data/try2_%d.csv'%sec, 'wb') as csvfile:
    #writer = csv.writer(csvfile)
  '''
  sky_data = SkyCoord(coo1, unit='deg', frame=FK5, equinox='J2000.0')
  gal = sky_data.transform_to(Galactic)
  coo1 = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)

  sky_data = SkyCoord(coo2, unit='deg', frame=FK5, equinox='J2000.0')
  gal = sky_data.transform_to(Galactic)
  coo2 = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
  '''

  if len2>len1:
    for i in range(len1):
      #print(i)
      tmp_co = coo2-coo1[i,:]
      mask = (np.absolute(tmp_co[:,0])<=bound) & (np.absolute(tmp_co[:,1])<=bound)
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
      if flux is not None:
        f = flux[i]/50.
      else:
        f = 1.
      weights = np.concatenate((weights,weight[mask]/f), axis=0)
  else:
    for i in range(len2):
      #print(i)
      tmp_co = coo2[i,:]-coo1
      mask = (np.absolute(tmp_co[:,0])<=bound) & (np.absolute(tmp_co[:,1])<=bound)
      tmp_co = tmp_co[np.absolute(tmp_co[:,0])<=bound,:]
      tmp_co = tmp_co[np.absolute(tmp_co[:,1])<=bound,:]
      co_rel = np.concatenate((co_rel, tmp_co), axis = 0)
      if flux is not None:
        f = flux[mask]/50.
      else:
        f = np.ones(tmp_co.shape[0])
      weights = np.concatenate((weights,np.repeat(weight[i], tmp_co.shape[0])/f), axis=0)

  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel[1:],1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                             bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]), weights=weights)
  print H.shape

  t = np.percentile(H, 90)
  t = 1
  print 'threshold:{0}'.format(t) 
  hp = H.copy()

  #plt.colorbar()
  '''
  plt.plot(co_rel[:,0],co_rel[:,1], '.k')
  plt.xlim(-0.02,0.02)
  plt.ylim(-0.02,0.02)
  '''
  #plt.savefig('/home/dw1519/galex/plots/co3/{0}.png'.format(idx))
  data = H.byteswap(True).newbyteorder()
  data = data.copy(order='C')
  data = data.byteswap(True).newbyteorder()
  c_data = c3.find_centroid(data, t)
  #print c_data
  if c_data is not None:
    #cy, cx, max_value, flux, x, y = c_data
    cy, cx, max_value, flux = c_data
    cy+=0.5
    cx+=0.5
    centroid = wcs.wcs_pix2world(wcs.sip_foc2pix([[cx, cy]],1),1)[0]
    #c_sex = wcs.wcs_pix2world(wcs.sip_foc2pix([[x+0.5, y+0.5]],1),1)[0]
    #print c_sex
  else:
    centroid = [0.,0.]
    max_value = 0
    flux = 500
  if centroid[0]>1:
    centroid[0] = centroid[0]-360.
  #if c_sex[0]>1:
  #  c_sex[0] = c_sex[0]-360.
  #print 'max:{0}'.format(max_value)
  #print cx,cy

  plt.imshow(np.log2(hp+10**-10),interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[bound, -bound, -bound, bound], origin='lower', vmin=0)
  #plt.imshow(hp+10**-10,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[0.01*3600, -0.01*3600, 0.01*3600, -0.01*3600], origin='upper', norm=LogNorm(10**0,np.max(hp)))  
  plt.colorbar()
  plt.xlim(-bound, bound)
  plt.ylim(-bound, bound)
  plt.ylabel('gb')
  plt.xlabel('gl')
  plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
  plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
  plt.plot(centroid[0],centroid[1],'+r', markersize=8)
  #plt.plot(c_sex[0],c_sex[1],'+g', markersize=8)
  plt.show()
  #plt.savefig('/home/dw1519/dw1519/galex/plots/co239_fix/{0}.pdf'.format(idx))
  #plt.clf()


  return centroid, max_value, flux

def get_data_fast(tranges, scan_name):
  look_up = load_obj('../data/%s/split/lookup'%scan_name)
  data = []
  trange = [0 , 0]
  num_interval = len(tranges)
  tranges = np.array(tranges, dtype='double')
  tranges = (tranges*1000).astype(int)

  initial_time = trange[0] = int(tranges[0,0]) #int(reader.next()[0]) + 0
  final_time = trange[1] = int(tranges[num_interval-1, 1])
  for i in range(1000):
    if '%d'%(initial_time) in look_up:
      print 'initial found:%d'%initial_time
      break
    initial_time += 5

  for i in range(1000):
    if '%d'%(final_time) in look_up:
      print 'final found:%d'%final_time
      break
    final_time -= 5

  for i in range(num_interval):
    data.append(np.array([], dtype='float64').reshape(0,11))

  if initial_time>=final_time:
    return data
  
  start_index = look_up['%d'%initial_time]
  end_index = look_up['%d'%final_time]
  csv_list = []
  for i in range(start_index, end_index+1):
    csv = np.loadtxt('../data/%s/split/%d.csv'%(scan_name, i),delimiter=',')
    time = csv[:,0]
    for j in range(num_interval):
      mask = (time>=tranges[j,0]) & (time<tranges[j,1])
      if np.sum(mask>0):
        data[j] = np.concatenate([data[j], csv[mask]], axis=0)
  return data

def angle_filter(data, skypos, angle):
  center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
  rad = np.cos(angle*np.pi/180.)
  #print(center, rad)

  data_rad = data*np.pi/180.
  X, Y, Z = pycoo.spherical_to_cartesian(1, data_rad[:,1], data_rad[:,0])
  data_car = np.array([X,Y,Z], dtype='float64').T

  sep = np.dot(data_car, center)
  mask = sep>=rad
  coo = data[sep>=rad, :]

  return coo, mask

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

def load_catalog():
    hdulist = pyfits.open('../data/tycho2.fits')
    star_data = hdulist[1].data
    length = star_data.shape[0]
    stars = np.zeros((length,2))
    stars[:,0] = star_data['RAJ2000']
    stars[:,1] = star_data['DEJ2000']
    stars_rad = stars*np.pi/180.
    X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
    stars_car = np.array([X,Y,Z], dtype='float64').T
    hdulist.close()

    return stars, stars_car, star_data

def load_catalog_matched(scan, cut=0):
    cfile='/scratch/dw1519/galex/data/extraction/Starcat-07-22-2017/starcat_{0}mapweight_fec_fwhm.npy'.format(scan)
    stars_data = np.load(cfile)
    stars = np.array([stars_data['RAJ2000'], stars_data['DEJ2000']]).T
    stars_rad = stars*np.pi/180.
    X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
    stars_car = np.array([X,Y,Z], dtype='float64').T
    flux = np.array([stars_data['flux'], stars_data['nuv']]).T
    cut = flux[:,1]>=cut

    return stars[cut], stars_car[cut], flux[cut]


def run_one_r_sec(pid, scan_name, step, resolution, duration, bound, asp_cal, start, end, dis_map, slope, inter, cut, tmp_dir, return_dict):

    print('run one r sec')

    num_co = int(resolution/step)
    
    #load asp file
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%scan_name)
    guide_time = hdulist[1].data['T']
    length = co_data.shape[0]
    hdulist.close()

    #stars, stars_car, star_data = load_catalog()
    scan = '{0:04d}'.format(int(re.split('_',scan_name)[3]))
    stars, stars_car, star_data = load_catalog_matched(scan, 0)#14.4)
    print stars.shape

    angle_list = [0.]

    print length
    start = 250#1300#493#1200
    #end = 1340
    skypos = [0.0, 0.0]
    skyrange = [2*bound, 2*bound]#[0.04, 0.04]
    for initial_sec in range(start, end):
      print 'time:%d'%initial_sec
      initial_time = guide_time[initial_sec]
      #rot = intitial_asp[3]
      #center = np.array([intitial_asp[1], intitial_asp[2]])
      #use the first order calibrated asp instead of the original one
      #center = cal_asp_r(offsets, cal_initial, intitial_asp[0]*1000, intitial_asp[1:3], 200) 
      #print(intitial_asp)

      tranges = []

      for sec in range(num_co):
        tranges.append([initial_time+step*sec-duration/2., initial_time+step*sec+duration/2.])
      print tranges
      time_c = np.mean(tranges, axis=1)
      print 'center time:'
      print time_c

      #data = get_data_cal(tranges, scan_name, cal_initial, offsets)
      data = get_data_fast(tranges, scan_name)
      print 'number of the photons:'
      print len(data)

      centroids = []
      center_time = []
      for sec in range(num_co):
        center_time.append(time_c[sec])
        if len(data[sec]) < duration*20000:
          centroid = np.array([0.0, 0.0, 0.0])
          centroids.append(centroid)
          print centroid
          continue
        arg = np.argmin(np.absolute(asp_cal[:,0]-time_c[sec]))
        center = asp_cal[arg, 1:3]

        data_sec, time_sec = dis_correct(np.array(data[sec], dtype='float64'), asp_cal, dis_map, slope, inter, cut)

        coo1 = np.array(data_sec, dtype='float64')
        aperture = 0.69
        #coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]
        #coo2 = catalog_fits.get_catalog_tycho(center, aperture)
        #coo2 = catalog_fits.get_catalog_tycho_fast(center, aperture, stars, stars_car, star_data)
        coo2, flux = catalog_fits.get_catalog_matched(center, aperture, stars, stars_car, star_data)

        coo1, mask = angle_filter(coo1, center, 0.6)
        time_sec = time_sec[mask]/1000.#np.array(data[sec], dtype='float64')[mask,0]/1000.
        print time_sec.shape
        print np.std(time_sec-time_c[sec])
        weights = gaussian(time_sec-time_c[sec], 0, 4*np.std(time_sec-time_c[sec]))

        #print coo2
        coo2, mask = angle_filter(coo2, center, 0.6)
        flux = flux[mask]

        max_now = -1000
        cent_now = []
        for angle in angle_list:
          centroid, max_value, flux = get_corr_map(coo2, coo1, skypos, skyrange, sec, 0.0002, initial_sec, '1', weights, flux[:,0])
          if max_value>max_now:
            max_now = max_value
            cent_now = np.append(centroid, angle) #centroid.append(angle)

        centroids.append(cent_now)
      print centroids
      np.save(tmp_dir+'/time_rot%d.npy'%(initial_sec), center_time)
      np.save(tmp_dir+'/centroids_rot%d.npy'%(initial_sec), centroids)


if __name__ == '__main__':

#hpc single, identify outliers, 0.5s, secondary co
  if True:

    parser = argparse.ArgumentParser(description='galex scan cross-correlation')
    parser.add_argument('scan_name', nargs=1, help="scan name")
    parser.add_argument('np', nargs=1, type=int, help="number of process")
    parser.add_argument('step', nargs=1, type=float, help="time spacing between two attitude solution in seconds")
    parser.add_argument('resolution', nargs=1, type=float, help="input asprta resolution")
    parser.add_argument('duration', nargs=1, type=float, help="length of photons to use in seconds")
    parser.add_argument('bound', nargs=1, type=float, help="bound of the centroids")
    parser.add_argument('asprta', nargs=1, help="path to the input asprta file")
    parser.add_argument('tmp_dir', nargs=1, help="directory to temperay files")
    parser.add_argument('out_suffix', nargs=1, help="suffix of the output file")
    parser.add_argument('-c', '--cut', action='store_true', help="whether to cut the photons by ya and q")
    parser.add_argument('-d', '--dis', nargs=2, help="path to the two distortion maps")
    parser.add_argument('-y', '--ya', nargs=1, help="path to the ya correction model")

    args = parser.parse_args()

    name = args.scan_name[0]

    asprta = args.asprta[0]
    
    p_num = args.np[0]

    step = args.step[0]
    resolution = args.resolution[0]
    duration = args.duration[0]
    bound = args.bound[0]

    tmp_dir = args.tmp_dir[0]
    out_suffix = args.out_suffix[0]

    cut = args.cut

    if args.dis is not None:
      dis_l = args.dis[0]
      dis_h = args.dis[1]
    else:
      dis_l = None
      dis_h = None

    if args.ya is not None:
      ya = args.ya[0]
    else:
      ya = None

    print("scan name: {0}".format(name))
    print("number of process: {0}".format(p_num))
    print("step: {0}".format(step))
    print("resolution: {0}".format(resolution))
    print("duration: {0}".format(duration))
    print("bound: {0}".format(bound))
    print("asprta: {0}".format(asprta))
    print("tmp dir: {0}".format(tmp_dir))
    print("output suffix: {0}".format(out_suffix))
    print("cut: {0}".format(cut))
    print("dis map low xa: {0}".format(dis_l))
    print("dis map high xa: {0}".format(dis_h))
    print("ya model: {0}".format(ya))

    hdulist = pyfits.open(asprta)
    co_data = hdulist[1].data
    #try:
    #  asp_cal = np.load('../data/photon_list/%s_asp_new.npy'%name)
    #except IOError:
    n=resolution/0.005 #100
    T = co_data['T']
    ra = co_data['ra']
    dec = co_data['dec']
    roll = co_data['roll']
    t_new = np.arange((T.shape[0]-1)*n)*0.005+T[0]
    ra_new = np.interp(t_new, T, ra)
    dec_new = np.interp(t_new, T, dec)
    roll_new = np.interp(t_new, T, roll)
    asp_cal = np.array([t_new, ra_new, dec_new, roll_new]).T
    #np.save('../data/photon_list/%s_asp_fix.npy'%name, asp_cal)
    hdulist.close()

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    guide_data = hdulist[1].data
    length = guide_data.shape[0]
    chunk_len = int(length/p_num)
    hdulist.close()

    test = tmp_dir+'/test'
    path = os.path.dirname(tmp_dir)
    if not os.path.exists(path):
        os.makedirs(path)

    if dis_l is not None:
      detector_size = [50,50]

      dis_map_l = np.load(dis_l)
      dis_map_h = np.load(dis_h)
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
    else:
      dis_map = None

    if ya is not None:
      ya_info = np.load(ya)
      slope = ya_info[0]
      inter = ya_info[1]
    else:
      slope = None
      inter = None

    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_r_sec, args=(pid, name, step, resolution, duration, bound, asp_cal, start, end, dis_map, slope, inter, cut, tmp_dir, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    #run_one_r_sec(pid, name, 1., start, end, return_dict)

    p = Process(target=run_one_r_sec, args=(pid, name, step, resolution, duration, bound, asp_cal, start, end, dis_map, slope, inter, cut, tmp_dir, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'

    #centroid_csv.generate_first_offsets(name)
    centroid_csv.generate_new_offsets_new(name, asprta, out_suffix, tmp_dir)

