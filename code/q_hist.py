import numpy as np
import matplotlib.pyplot as plt
import sys
from fast_histogram import histogram2d
from scipy import interpolate
import gnomonic as gn
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
import glob
import re

def dis_correct(data, dis_map=None, ya_corr=None, q_corr=None, cut=False):
    if cut:
      ya_mask = data[:,2]>=2
      q_mask = data[:,3]>=0
      data = data[ya_mask&q_mask]

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

    #ra, dec = gn.gnomrev_simple(data[:,4], data[:,5], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return data[:,4:6]

def dis_correct_c(data, dis_map=None, ya_corr=None, q_corr=None, cut=False):
    if cut:
      ya_mask = data[:,2]>=2
      q_mask = data[:,3]>5
      data = data[ya_mask&q_mask]

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

    #ra, dec = gn.gnomrev_simple(data[:,4], data[:,5], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return data[:,4:6]



def detector2gon(pos_list, size):
  pos_array = np.array(pos_list, dtype='float64')
  pos = (pos_array/50.-0.5)*36000*800*0.001666
  return pos


if __name__ == '__main__':
    scan = sys.argv[1]
    dis_l = sys.argv[2]
    dis_h  = sys.argv[3]
    ya = sys.argv[4]
    q = sys.argv[5]
    data_path = sys.argv[6]

    num= int(re.split('_',scan)[3])
    print num
    if num>=5 and num<=185:
        dis_l = '/scratch/dw1519/galex/data/distortion/5-185-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/5-185-xa_high-centroids.npy'
    elif num>=194 and num<=437:
        dis_l = '/scratch/dw1519/galex/data/distortion/194-437-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/194-437-xa_high-centroids.npy'        
    elif num>=446 and num<=887:
        dis_l = '/scratch/dw1519/galex/data/distortion/446-887-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/446-887-xa_high-centroids.npy'        
        #elif num>=896 and num<=1157:
    #    content = re.sub(r'^dis=.*$', 'dis=/scratch/dw1519/galex/data/distortion/896-1157', content, flags=re.MULTILINE)
    elif num>=1166 and num<=1373:
        dis_l = '/scratch/dw1519/galex/data/distortion/1166-1373-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/1166-1373-xa_high-centroids.npy'        
    else:
        dis_l = '/scratch/dw1519/galex/data/distortion/194-437-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/194-437-xa_high-centroids.npy' 

    pixsz = 0.0005555555556
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

    if q is not None:
        q_corr = np.load(q+'/'+scan+'/q_corr.npy')
    else:
        q_corr = None
    if ya is not None:
        ya_corr = np.load(ya+'/'+scan+'/ya_corr.npy')
    else:
        ya_corr = None

    csvs = glob.glob("{0}/{1}/split/*.csv".format(data_path, scan))

    size = 800
    detector = np.zeros((size,size))
    detector_c = np.zeros((size,size))
    cut = True
    for csv_file in csvs:
        print 'loding {0}'.format(csv_file)
        #sdata = np.loadtxt(csv_file, delimiter=',', usecols=[0,8,9])
        data = np.loadtxt(csv_file, delimiter=',', usecols=[0,3,4,5,6,7])
        print 'loaded,shape:{0}'.format(data.shape)

        data_f = dis_correct(data, dis_map, ya_corr, q_corr, cut)
        pos = ((data_f/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*size
        #print np.max(pos)
        #print np.min(pos)
        H=histogram2d(pos[:,0], pos[:,1],\
                                  bins=[size,size], range=([ [0,size],[0,size] ]))
        detector += H

        data_c = dis_correct_c(data, dis_map, ya_corr, q_corr, cut)
        pos = ((data_c/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*size
        #print np.max(pos)
        #print np.min(pos)
        H=histogram2d(pos[:,0], pos[:,1],\
                                  bins=[size,size], range=([ [0,size],[0,size] ]))
        detector_c += H

    np.save("{0}/{1}/full.npy".format(data_path, scan), detector)
    np.save("{0}/{1}/cut.npy".format(data_path, scan), detector_c)

