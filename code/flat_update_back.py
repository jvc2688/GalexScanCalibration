import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import sys
#from fast_histogram import histogram2d
from scipy import interpolate
import gnomonic as gn
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
import glob
import load_data
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits as pyfits
from astropy.nddata.utils import block_reduce
import scipy

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
    return np.array([data[:,0],ra,dec, data[:,4], data[:,5]]).T

def detector2gon(pos_list, size):
  pos_array = np.array(pos_list, dtype='float64')
  pos = (pos_array/50.-0.5)*36000*800*0.001666
  return pos


if __name__ == '__main__':
    size=800
    resample = 2
    outdir = 'back_iterate1'
    f_scans = ['AIS_GAL_SCAN_00005_0001', 'AIS_GAL_SCAN_00446_0001', 'AIS_GAL_SCAN_00995_0001',
               'AIS_GAL_SCAN_01247_0001', 'AIS_GAL_SCAN_02120_0002',
                'AIS_GAL_SCAN_03326_0002']
    full = np.zeros((size,size))
    cut = np.zeros((size,size))
    for scan in f_scans:
        full += np.load('/beegfs/dw1519/galex/data/'+scan+'/full.npy')
        cut += np.load('/beegfs/dw1519/galex/data/'+scan+'/cut.npy')
    flat = pyfits.open('/scratch/dw1519/galex/data/cal/NUV_flat.fits')[0].data
    #factor = 4
    #flat = block_reduce(flat, factor, np.mean)
    frac = cut/full
    frac[np.isnan(frac)] = 0.
    flat = frac*flat
    print flat.shape
    np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat_in.npy', flat)
    plt.imshow(flat, vmin=0, vmax=1.5)
    plt.colorbar()
    plt.savefig('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat_in.pdf', dpi=190)
    plt.clf()
    scans = ['0023', '0212', '0302', '0464', '0500', '0689', '0806', '0941', '1247', '1508', '1778', 
            '2120', '2291', '2453','2759', '2930', '3155', '3227', '3326', '3497']
    date = '08-21-2017'#'08-21-2017'
    data = np.zeros((size,size))
    flat = np.load('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat0.npy')
    data = np.load('/scratch/dw1519/galex/data/star_photon/'+outdir+'/data.npy')
    data = block_reduce(data, resample)
    N=10
    for i in range(1, N):
        model = np.zeros((size,size))
        A_list = []
        for scan in scans:
            with open('../name_scan/%s'%scan) as f:
                name_list = f.read().splitlines()
            print name_list
            if i==0:
                data += np.load('/scratch/dw1519/galex/data/star_photon/back/'+scan+'_data.npy')
            exp = np.load('/scratch/dw1519/galex/data/star_photon/back/'+scan+'_exp.npy')
            if len(exp.shape)<2:
                exp = np.concatenate(exp,axis=0)
            print exp.shape
            star_list = np.load('/scratch/dw1519/galex/data/star_photon/back/'+scan+'_bstar.npy')
            print star_list.shape
            stars, ct = np.unique(exp[:,0], return_counts=True)
            print stars.shape
            star_mask = (ct>10)
            print np.sum(star_mask)
            stars = stars[star_mask]
            for star in stars:
                print star
                star_mask = exp[:,0]==star
                p = exp[star_mask]
                x = p[:,1].astype(int)
                y = p[:,2].astype(int)
                A = np.sum(p[:,-2])/np.dot(p[:,-3],flat[x,y])
                A_list.append(A)
                model[x,y] += A*p[:,-3]        
        if i==0:
            np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/data.npy', data)
            data = block_reduce(data, resample)
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/model{0}.npy'.format(i), model)
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/A{0}.npy'.format(i), A_list)
        model = block_reduce(model, resample)
        model[model==0] = 1.
        flat = data/model
        flat = scipy.ndimage.interpolation.zoom(flat, resample, order=0,prefilter=False)
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat{0}.npy'.format(i), flat)
        plt.imshow(flat, vmin=0, vmax=1.5)
        plt.colorbar()
        plt.savefig('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat{0}.pdf'.format(i), dpi=190)
        plt.clf()
