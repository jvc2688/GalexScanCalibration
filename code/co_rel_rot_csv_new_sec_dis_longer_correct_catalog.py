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
import centroid_plot_csv_new_2 as centroid_csv
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
from astropy import units as u
import load_data
import glob


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def dis_correct(data, asp, dis_map=None, ya_corr=None, q_corr=None, cut=False):
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

    if q_corr is not None:
      q_mask = data[:,5]<24
      data[q_mask,6] -= q_corr[data[q_mask,5].astype(int),-2]*36000*800*0.001666/2400.
      data[q_mask,7] -= q_corr[data[q_mask,5].astype(int),-1]*36000*800*0.001666/2400.

    if ya_corr is not None:
      data[:,6] -= ya_corr[data[:,4].astype(int)-2,-2]*36000*800*0.001666/2400.
      data[:,7] -= ya_corr[data[:,4].astype(int)-2,-1]*36000*800*0.001666/2400.

    ra, dec = gn.gnomrev_simple(data[:,6], data[:,7], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return np.array([ra,dec]).T, data[:,0]


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

def get_corr_map(pid, coo1, coo2, skypos, skyrange, pixsz, time, suffix, weight, flux=None):
  imsz = imagetools.deg2pix(skypos, skyrange, pixsz).astype(int)
  count = np.zeros(imsz)
  bound = skyrange[0]/2.

  wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)

  catalog = SkyCoord(coo1, unit='deg', frame='fk5', equinox='J2000.0')
  sky_data = SkyCoord(coo2, unit='deg', frame='fk5', equinox='J2000.0')
  print '{1}:catalog:{0}'.format(catalog.shape[0], pid)
  print '{1}:data:{0}'.format(sky_data.shape[0], pid)
  idxc, idxcatalog, d2d_p, d3d_p = sky_data.search_around_sky(catalog, 0.0045*u.deg)
  print '{1}:match:{0}'.format(idxc.shape[0], pid)

  if idxc.shape[0]<10:
    centroid = [0.,0.]
    max_value = 0
    flux = 500
    return centroid, max_value, flux, idxc.shape[0]

  co_rel = coo2[idxcatalog] - coo1[idxc]
  if flux is not None:
    weights = weight[idxcatalog]/flux[idxc]*50.
  else:
    weights = weight[idxcatalog]

  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co_rel,1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                             bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]), weights=weights)

  t = np.percentile(H, 90)
  t = 1
  print '{1}:threshold:{0}'.format(t, pid) 
  #hp = H.copy()

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
  '''
  plt.imshow(np.log2(hp+10**-10),interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[bound, -bound, -bound, bound], origin='lower', vmin=0)
  #plt.imshow(hp+10**-10,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[0.01*3600, -0.01*3600, 0.01*3600, -0.01*3600], origin='upper', norm=LogNorm(10**0,np.max(hp)))  
  plt.colorbar()
  plt.xlim(-bound/3, bound/3)
  plt.ylim(-bound/3, bound/3)
  plt.ylabel('gb')
  plt.xlabel('gl')
  plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
  plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
  plt.plot(centroid[0],centroid[1],'+r', markersize=8)
  #plt.plot(c_sex[0],c_sex[1],'+g', markersize=8)
  plt.show()
  #plt.savefig('/home/dw1519/dw1519/galex/plots/co239_fix/{0}.pdf'.format(idx))
  #plt.clf()
  '''
  return centroid, max_value, flux, idxc.shape[0]

def get_data_fast(tranges, scan_name):
  look_up = load_obj('/beegfs/dw1519/galex/data/%s/split/lookup'%scan_name)
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
    csv = np.loadtxt('/beegfs/dw1519/galex/data/%s/split/%d.csv'%(scan_name, i),delimiter=',')
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

def convert_catalog(cfile):
    tycho2 = pyfits.open('../data/tycho2.fits')[1].data

    sfile=re.split('\.',cfile)[0]+'.txt'
    print(sfile)
    try:
      df = load_data.load_catalog(sfile)
    except IOError:
      print('skip')
      return None
    c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
    catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask=d2d<0.001*u.degree
    print(np.sum(mask))
    dtype = np.dtype([('tycho_num', int), ('Glon', '>f4'), ('Glat', '>f4'), ('RAJ2000', '>f4'), ('DEJ2000', '>f4'),
             ('flux', float), ('nuv', float), ('gl', float), ('gb', float)])
    matched_tycho = tycho2[idx[mask]]
    matched_df = df[mask]
    matched_catalog = np.core.records.fromarrays(np.array([idx[mask], matched_tycho['Glon'], matched_tycho['Glat'],
            matched_tycho['RAJ2000'], matched_tycho['DEJ2000'], np.array(matched_df['FLUX_AUTO']),
            np.array(matched_df['nuv']), np.array(matched_df['gl']),
            np.array(matched_df['gb'])]), dtype=dtype)
    print(matched_catalog.shape)
    np.save(cfile, matched_catalog)
    return matched_catalog

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
    cfile='/scratch/dw1519/galex/fits/scan_map/catalog/08-21-2017_1/starcat_{0}mapweight_fec_fwhm.npy'.format(scan)
    try:
      stars_data = np.load(cfile)
    except IOError:
      stars_data = convert_catalog(cfile)
    stars = np.array([stars_data['RAJ2000'], stars_data['DEJ2000']]).T
    stars_rad = stars*np.pi/180.
    X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
    stars_car = np.array([X,Y,Z], dtype='float64').T
    flux = np.array([stars_data['flux'], stars_data['nuv']]).T
    cut = flux[:,1]>=cut

    return stars[cut], stars_car[cut], flux[cut]


def run_one_r_sec(pid, scan_name, step, resolution, duration, bound, asp_cal, time_list, dis_map, ya_corr, q_corr, cut, tmp_dir, return_dict):

    print('run one r sec')

    num_co = int(resolution/step)
    
    #stars, stars_car, star_data = load_catalog()
    scan = '{0:04d}'.format(int(re.split('_',scan_name)[3]))
    stars, stars_car, star_data = load_catalog_matched(scan, 14.4)
    print stars.shape

    angle_list = [0.]


    skypos = [0.0, 0.0]
    skyrange = [2*bound, 2*bound]
    centroids = []
    n_match = []
    for time_c in time_list:
      tranges = [[time_c-duration/2., time_c+duration/2.]]
      #print tranges
      print '{1}:center t:{0}'.format(time_c, pid)

      #data = get_data_cal(tranges, scan_name, cal_initial, offsets)
      data_sec = get_data_fast(tranges, scan_name)[0]
      #print 'number of the photons:'
      #print len(data_sec)

      if len(data_sec) < duration*15000:
        centroid = np.array([0.0, 0.0, 0.0])
        print '{1}:cent:{0}'.format(centroid, pid)
        centroids.append(centroid)
        nm=0
        n_match.append(nm)
        continue
      arg = np.argmin(np.absolute(asp_cal[:,0]-time_c))
      center = asp_cal[arg, 1:3]

      data_sec, time_sec = dis_correct(np.array(data_sec, dtype='float64'), asp_cal, dis_map, ya_corr, q_corr, cut)

      coo1 = np.array(data_sec, dtype='float64')
      aperture = 0.69

      coo2, flux = catalog_fits.get_catalog_matched(center, aperture, stars, stars_car, star_data)

      coo1, mask = angle_filter(coo1, center, 0.6)
      time_sec = time_sec[mask]/1000.#np.array(data[sec], dtype='float64')[mask,0]/1000.
      weights = gaussian(time_sec-time_c, 0, step)

      #print coo2
      coo2, mask = angle_filter(coo2, center, 0.6)
      flux = flux[mask]
      if len(coo2)<=0:
        centroid = np.array([0.0, 0.0, 0.0])
        print '{1}:cent:{0}'.format(centroid, pid)
        centroids.append(centroid)
        nm=0
        n_match.append(nm)
        continue
      max_now = -1000
      cent_now = []
      for angle in angle_list:
        centroid, max_value, flux, nm = get_corr_map(pid, coo2, coo1, skypos, skyrange, 0.0002, time_c, '1', weights, flux[:,0])
        if max_value>max_now:
          max_now = max_value
          cent_now = np.append(centroid, angle) #centroid.append(angle)
      print '{1}:cent:{0}'.format(cent_now, pid)
      centroids.append(cent_now)
      n_match.append(nm)
    #print centroids
    np.save(tmp_dir+'/centroids_tmp%d.npy'%(pid), centroids)
    np.save(tmp_dir+'/match_tmp%d.npy'%(pid), n_match)


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
    parser.add_argument('-q', '--q', nargs=1, help="path to the q correction model")

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

    if args.q is not None:
      q = args.q[0]
    else:
      q = None

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
    print("q model: {0}".format(q))

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

    '''
    ra_new = np.interp(t_new, T, ra)
    dec_new = np.interp(t_new, T, dec)
    roll_new = np.interp(t_new, T, roll)
    '''
    spl_ra = splrep(T, ra)
    spl_dec = splrep(T, dec)
    spl_roll = splrep(T, roll)
    ra_new = splev(t_new, spl_ra)
    dec_new = splev(t_new, spl_dec)
    roll_new = splev(t_new, spl_roll)
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
      ya_corr = np.load(ya)
    else:
      ya_corr = None

    if q is not None:
      q_corr = np.load(q)
    else:
      q_corr = None

    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    guide_data = hdulist[1].data
    guide_time = guide_data['T']
    hdulist.close()
    guide_res = 1.

    scst_file = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%name)
    scst_data = scst_file[1].data
    scst_time = scst_data['pktime'].copy()
    hv = scst_data['hvnom_nuv'].copy()
    scst_file.close()

    nstep = int(guide_res/step)
    time_c = np.repeat(guide_time, nstep)
    for n in range(nstep):
      time_c[n::nstep] += n*step

    scst_ix = np.digitize(time_c, scst_time)-1
    ix_mask = scst_ix<scst_time.shape[0]
    scst_ix = scst_ix[ix_mask]
    time_c = time_c[ix_mask]
    time2run = time_c[hv[scst_ix]>0]

    length = time2run.shape[0]
    chunk_len = int(length/p_num)


    manager = Manager()
    return_dict = manager.dict()
    p_list = []
    pid = 0
    for pid in range(p_num-1):
      start = pid*chunk_len
      end = (pid+1)*chunk_len
      p = Process(target=run_one_r_sec, args=(pid, name, step, resolution, duration, bound, asp_cal, time2run[start:end], dis_map, ya_corr, q_corr, cut, tmp_dir, return_dict))
      p.start()
      p_list.append(p)

    pid = p_num-1
    start = pid*chunk_len
    end = length
    #run_one_r_sec(pid, name, 1., start, end, return_dict)

    p = Process(target=run_one_r_sec, args=(pid, name, step, resolution, duration, bound, asp_cal, time2run[start:end], dis_map, ya_corr, q_corr, cut, tmp_dir, return_dict))
    p.start()
    p_list.append(p)

    for p in p_list:
      p.join()
    print 'all done'

    offset = np.zeros((time_c.shape[0],3))
    n_match = np.zeros(time_c.shape[0])
    centers = []
    nms = []
    for i in range(p_num):
        centers.append(np.load(tmp_dir+'/centroids_tmp%d.npy'%(i)))
        nms.append(np.load(tmp_dir+'/match_tmp%d.npy'%(i)))
    centroids = np.concatenate(centers, axis=0)
    nms = np.concatenate(nms, axis=0)
    offset[hv[scst_ix]>0] = centroids
    n_match[hv[scst_ix]>0] = nms

    np.save(tmp_dir+'/offsets_%s.npy'%(out_suffix), offset)
    np.save(tmp_dir+'/time_%s.npy'%(out_suffix), time_c)
    np.save(tmp_dir+'/match_%s.npy'%(out_suffix), n_match)

    tmp_files = glob.glob(tmp_dir+"/centroids_tmp*")
    for tmp_file in tmp_files:
        os.remove(tmp_file)

    tmp_files = glob.glob(tmp_dir+"/time_tmp*")
    for tmp_file in tmp_files:
        os.remove(tmp_file)

    tmp_files = glob.glob(tmp_dir+"/match_tmp*")
    for tmp_file in tmp_files:
        os.remove(tmp_file)

    hv_mask = hv[scst_ix]>0
    nan_mask = np.zeros(time_c.shape[0], dtype=bool)
    nan_mask[hv_mask] = (np.sum(centroids,axis=1)==0)
    #centroid_csv.generate_first_offsets(name)
    centroid_csv.generate_new_offsets_new(name, asprta, out_suffix, tmp_dir, p_num, hv_mask, nan_mask)
