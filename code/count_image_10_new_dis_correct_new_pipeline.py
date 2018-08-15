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
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
import argparse
from fast_histogram import histogram2d


def dis_correct(data, asp, dis_map=None, ya_corr=None, q_corr=None, cut=False):
    if cut:
      ya_mask = data[:,2]>=2
      q_mask = data[:,3]>5
      data = data[ya_mask&q_mask]
    ix = np.digitize(data[:,0]/1000., asp[:,0])-1

    if dis_map is not None:
      low_mask = data[:,1]<16
      high_mask = data[:,1]>=16

      coo = data[low_mask, 4:6]
      data[low_mask, 4] -= dis_map['xi_l'](coo)
      data[low_mask, 5] -= dis_map['eta_l'](coo)

      coo = data[high_mask, 4:6]
      data[high_mask, 4] -= dis_map['xi_h'](coo)
      data[high_mask, 5] -= dis_map['eta_h'](coo)

    if q_corr is not None:
      q_mask = data[:,3]<24
      data[q_mask,4] -= q_corr[data[q_mask,3].astype(int),-2]*36000*800*0.001666/2400.
      data[q_mask,5] -= q_corr[data[q_mask,3].astype(int),-1]*36000*800*0.001666/2400.

    if ya_corr is not None:
      data[:,4] -= ya_corr[data[:,2].astype(int)-2,-2]*36000*800*0.001666/2400.
      data[:,5] -= ya_corr[data[:,2].astype(int)-2,-1]*36000*800*0.001666/2400.

    ra, dec = gn.gnomrev_simple(data[:,4], data[:,5], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
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

  while total_p_num>0:
    csv_info_list = np.sort(csv_info_list, order='num_p')[::-1]
    csv_info_list[0]['p_num'] += 1
    csv_info_list[0]['num_p'] = float(csv_info_list[0]['num'])/csv_info_list[0]['p_num']
    total_p_num -= 1
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
  foc = wcs.all_world2pix(data,1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

def hist_fast(data, wcs, imsz):
  foc = wcs.all_world2pix(data,1)
  H=histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

def run_one_inter_fast(pid, file_name, imsz, wcs, csv_list, out_path, asp, scst_time, hv, dis_map, ya_corr, q_corr, cut, return_dict):
  count = np.zeros(imsz)
  step = 200
  tranges = []
  trange = [0, 0]

  i=0
  for csv_file in csv_list:
    print '{0}:loding {1}'.format(pid,csv_file)
    #sdata = np.loadtxt(csv_file, delimiter=',', usecols=[0,8,9])
    data = np.loadtxt(csv_file, delimiter=',', usecols=[0,8,9])
    print '{0}:loaded,shape:{1}'.format(pid,data.shape)
    if len(data.shape)<2:
      print '{0}:skip'.format(pid)
      continue
    #trim the data with scst flag
    scst_ix = np.digitize(data[:,0]/1000., scst_time)
    ix_mask = scst_ix<scst_time.shape[0]
    scst_ix = scst_ix[ix_mask]
    data = data[ix_mask]
    data = data[hv[scst_ix]>0]
    print '{0}:filter,shape:{1}'.format(pid,data.shape)
    if len(data)==0:
      print '{0}:skip'.format(pid)
      continue

    if i==0:
      trange[0] = float(data[0,0])/1000.
    trange[1] = float(data[-1,0])/1000.

    #data = dis_correct(data, asp, dis_map, ya_corr, q_corr, cut)
    print '{0}:correct,shape:{1}'.format(pid,data.shape)
    print '{0}:convert2gal'.format(pid)
    sky_data = SkyCoord(data[:,-2:], unit='deg', frame=FK5, equinox='J2000.0')
    gal = sky_data.transform_to(Galactic)
    data = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
    
    print '{0}:fast hist'.format(pid) 
    H = hist_fast(data, wcs, imsz)
    print '{0}:add'.format(pid) 
    count += H

    sky_data=gal=data=H = None
    gc.collect()
    print '{0}:{1}'.format(pid,i)
    i+=1
    if (i%50)==0:
      with open('{0}/{1}_gal_sec_count_tmp{2}.dat'.format(out_path, name_file, pid),'w') as f:
        f.write('%d'%i)
  tranges.append(trange)
  return_dict[pid] = (count, tranges)
  #np.save('{0}/{1}_gal_sec_count_tmp{2}.npy'.format(out_path, name_file, pid), count)
  #return_dict[pid] = tranges


if __name__ == '__main__':

#process interpolate fast
  if True:
    parser = argparse.ArgumentParser(description='galex scan count map')
    parser.add_argument('scan_name', nargs=1, help="scan name")
    parser.add_argument('np', nargs=1, type=int, help="number of process")
    parser.add_argument('data_path', nargs=1, help="path to the data")
    parser.add_argument('guide_path', nargs=1, help="path to the guide asp file")
    parser.add_argument('guide_suffix', nargs=1, help="suffix of the guide asp file")
    parser.add_argument('asprta', nargs=1, help="path to the asprta file")
    parser.add_argument('asprta_suffix', nargs=1, help="suffix the asprta file")
    parser.add_argument('asprta_resolution', nargs=1, type=float, help="resolution of asprta")
    parser.add_argument('scst', nargs=1, help="path to the scst file")
    parser.add_argument('out_path', nargs=1, help="path to the output directory")
    parser.add_argument('-c', '--cut', action='store_true', help="whether to cut the photons by ya and q")
    parser.add_argument('-d', '--dis', nargs=2, help="path to the two distortion maps")
    parser.add_argument('-y', '--ya', nargs=1, help="path to the ya correction model")
    parser.add_argument('-q', '--q', nargs=1, help="path to the q correction model")

    args = parser.parse_args()

    name_file = args.scan_name[0]

    total_p_num = args.np[0]

    data_path = args.data_path[0]

    guide_path = args.guide_path[0]
    guide_suffix = args.guide_suffix[0].lstrip()

    asprta = args.asprta[0]
    asprta_suffix = args.asprta_suffix[0].lstrip()
    resolution = args.asprta_resolution[0]
    scst = args.scst[0]
    
    out_path = args.out_path[0]

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

    print("scan name: {0}".format(name_file))
    print("number of process: {0}".format(total_p_num))
    print("guide asp path: {0}".format(guide_path))
    print("guide asp suffix: {0}".format(guide_suffix))
    print("asprta: {0}".format(asprta))
    print("asprta suffix: {0}".format(asprta_suffix))
    print("resolution: {0}".format(resolution))
    print("output: {0}".format(out_path))
    print("cut: {0}".format(cut))
    print("dis map low xa: {0}".format(dis_l))
    print("dis map high xa: {0}".format(dis_h))
    print("ya model: {0}".format(ya))
    print("q model: {0}".format(q))


    pixsz = 0.0005555555556 #0.000416666666666667 #0.00166667

    with open('../name_scan/%s'%name_file) as f:
      name_list = f.read().splitlines()
    print name_list

    csv_lists = {}
    csv_info_list = []
    q_corr_list = {}
    ya_corr_list = {}
    total_num_csv = 0
    gal_l = []
    dtype1 = np.dtype([('name', np.str_, 23), ('num', int), ('p_num', int), ('num_p', float)])
    for name in name_list:
      csvs = glob.glob("{0}/{1}/split/*.csv".format(data_path, name))
      csv_lists[name] = csvs
      num = len(csvs)
      csv_info_list.append((name, num, 1, num/1))
      total_num_csv += num

      if q is not None:
        q_corr_list[name] = np.load(q+'/'+name+'/q_corr.npy')
      else:
        q_corr_list[name] = None
      if ya is not None:
        ya_corr_list[name] = np.load(ya+'/'+name+'/ya_corr.npy')
      else:
        ya_corr_list[name] = None

      guide_asp = np.load('{0}/{1}_asp.npy'.format(guide_path, name))#np.load('../data/photon_list/%s_asp.npy'%name)
      sky_data = SkyCoord(guide_asp[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
      gal = sky_data.transform_to(Galactic)
      guide = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
      print np.min(guide[:,0])
      print np.max(guide[:,0])
      print asprta
      gal_l.append(np.mean(guide[:,0]))

    skypos = [np.mean(gal_l), 0.]
    print skypos
    print np.max(gal_l)
    print np.min(gal_l)
    skyrange = [1.55, 21.5]
    print skyrange
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz).astype(int)
    print imsz
    wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
    count = np.zeros(imsz)
    print skypos, skyrange

    #dis_map_l = np.load('../data/distortion/194-437-xa_low-centroids.npy')
    #dis_map_h = np.load('../data/distortion/194-437-xa_high-centroids.npy')
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

    '''
    if ya is not None:
      ya_corr = np.load(ya)
    else:
      ya_corr = None

    if q is not None:
      q_corr = np.load(q)
    else:
      q_corr = None
    '''

    csv_info_list = np.array(csv_info_list, dtype=dtype1)

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

      #try:
      #  asp = np.load('../data/photon_list/{0}{1}-correct_asp_new.npy'.format(name, suffix))
      #except IOError:
      co_data = pyfits.open('{0}/{1}{2}-asprta.fits'.format(asprta, name, asprta_suffix))[1].data
      n = resolution/0.005
      T = co_data['T']
      ra = co_data['ra']
      dec = co_data['dec']
      roll = co_data['roll']
      t_new = np.arange((T.shape[0]-1)*n)*0.005+T[0]
      '''
      f_ra = interp1d(T, ra, kind='cubic')
      f_dec = interp1d(T, dec, kind='cubic')
      f_roll = interp1d(T, roll, kind='cubic')

      ra_new = f_ra(t_new)
      dec_new = f_dec(t_new)
      roll_new = f_roll(t_new)
      '''

      spl_ra = splrep(T, ra)
      spl_dec = splrep(T, dec)
      spl_roll = splrep(T, roll)
      ra_new = splev(t_new, spl_ra)
      dec_new = splev(t_new, spl_dec)
      roll_new = splev(t_new, spl_roll)

      '''
      ra_new = np.interp(t_new, T, ra)
      dec_new = np.interp(t_new, T, dec)
      roll_new = np.interp(t_new, T, roll)
      '''
      asp = np.array([t_new, ra_new, dec_new, roll_new]).T
      #  np.save('../data/photon_list/{0}{1}-correct_asp_new.npy'.format(name, suffix), asp)

      scst_file = pyfits.open('{0}/{1}-scst.fits'.format(scst, name))
      scst_data = scst_file[1].data
      scst_time = scst_data['pktime'].copy()
      hv = scst_data['hvnom_nuv'].copy()
      scst_file.close()

      ya_corr = ya_corr_list[name]
      q_corr = q_corr_list[name]

      csvs = csv_lists[name]
      csvs = sorted(csvs, key=get_key)

      chunks = split_seq(csvs, p_num)
      #print chunks
      for chunk in chunks:
        print(len(chunk))
        p = Process(target=run_one_inter_fast, args=(pid, name, imsz, wcs, chunk, out_path,
                                                 asp, scst_time, hv, dis_map, ya_corr, q_corr, cut, return_dict))
        p.start()
        p_list.append(p)
        pid += 1
    print(pid)

    for p in p_list:
      p.join()
    print('all done')

    for i in range(0, pid):
      count += return_dict[i][0]#np.load('{0}/{1}_gal_sec_count_tmp{2}.npy'.format(out_path, name_file, i)) #return_dict[i][0]
      tranges += return_dict[i][1]

    print(tranges)
    print(count.shape)
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('{0}/count_map_{1}_count.fits'.format(out_path, name_file), clobber=False)

    for i in range(pid):
      #os.remove('{0}/{1}_gal_sec_count_tmp{2}.npy'.format(out_path, name_file, i))
      os.remove('{0}/{1}_gal_sec_count_tmp{2}.dat'.format(out_path, name_file, i))


