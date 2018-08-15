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
from multiprocessing import Process
from multiprocessing import Manager
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


def update(pid, stars, exp, size, flat, flatb, return_dict):
    data = np.zeros((size, size))
    model = np.zeros((size, size))
    A_list = []
    B_list = []
    nc = 5
    for star in stars:
        print star
        star_mask = exp[:,0]==star
        p = exp[star_mask]
        x = p[:,1].astype(int)
        y = p[:,2].astype(int)
        cts = p[:,-2]
        if np.sum(cts)<50:
          continue
        ERs = np.multiply(p[:,-3], flat[x,y])
        ERBs = np.multiply(p[:,-3], flatb[x,y])

        mask = ERs>0
        p = p[mask]
        length = p.shape[0]
        ERs = ERs[mask]
        ERBs = ERBs[mask]
        x = p[:,1].astype(int)
        y = p[:,2].astype(int)
        cts = p[:,-2]
        if np.sum(cts)<50:
          continue
        idx = np.random.randint(0, nc, length)
        C = np.bincount(idx, cts)
        while np.min(C)<5:
          idx = np.random.randint(0, nc, length)
          C = np.bincount(idx, weights=cts)
        ER = np.bincount(idx, weights=ERs)
        ERB = np.bincount(idx, weights=ERBs)
        w,r = scipy.optimize.nnls(np.vstack([ER,ERB]).T, C)
        #A, B = np.linalg.lstsq(np.vstack([ER,ERB]).T, C)[0]
        A = w[0]
        B = w[1]
        A_list.append(A)
        B_list.append(B)
        model[x,y] += A*p[:,-3]
        data[x,y] += (cts - B*ERBs)

    return_dict[pid] = (A_list, B_list, model, data)


if __name__ == '__main__':
    size=800
    outdir = 'iterate5'
    num_t = 5
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
    scans = ['0005', '0023', '0050', '0158', '0446', '0464', '0482',
             '0806', '0689','0941', '0995','1247','1508', '1778', 
             '3326', '2120', '2633','3047', '0212','0302','0500',
             '0743','2291','2453', '2759', '2885', '2930']#, '3155', '3227', '3497']
    date = '08-21-2017'#'08-21-2017'
    #flat = np.load('/scratch/dw1519/galex/data/star_photon/iterate/flat3.npy')
    flatb = np.load('/scratch/dw1519/galex/data/star_photon/back_iterate1/flat9.npy')
    N=10
    for i in range(N):
        data = np.zeros((size,size))
        model = np.zeros((size,size))
        A_list = []
        B_list = []
        for scan in scans:
            with open('/scratch/dw1519/galex/data/star_photon/'+outdir+'/now.dat', 'w') as f:
                f.write(scan)
            with open('../name_scan/%s'%scan) as f:
                name_list = f.read().splitlines()
            print name_list
            '''
            photon_list = []
            for name in name_list:
                photon_list.append(np.load('/scratch/dw1519/galex/data/star_photon/new/'+name+'_photon.npy'))
            photon_list = np.concatenate(photon_list, axis=0)
            '''
            exp = np.load('/scratch/dw1519/galex/data/star_photon/new/'+scan+'_exp.npy')
            if len(exp.shape)<2:
                exp = np.concatenate(exp,axis=0)
            print exp.shape
            star_list = np.load('/scratch/dw1519/galex/data/star_photon/nocut/'+name_list[0]+'_star.npy')
            print star_list.shape
            stars, ct = np.unique(exp[:,0], return_counts=True)#np.unique(photon_list[:,0], return_counts=True)
            nuv = star_list[stars.astype(int),-1]
            print stars.shape
            star_mask = (ct>10)&(nuv<17.5)&(nuv>15.5)#(ct>10) & (nuv<18) &(nuv>16.)
            print np.sum(star_mask)
            stars = stars[star_mask]

            stars_list = np.array_split(stars, num_t)
            manager = Manager()
            return_dict = manager.dict()
            p_list = []

            for pid in range(0, num_t):
                p = Process(target=update, args=(pid, stars_list[pid], exp, size, flat, flatb, return_dict))
                p.start()
                p_list.append(p)

            for p in p_list:
                p.join()

            for pid in range(0, num_t):
                result = return_dict[pid]
                A_list += result[0]
                B_list += result[1]
                model += result[2]
                data += result[3]

        #if i==0:
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/data{0}.npy'.format(i), data)
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/model{0}.npy'.format(i), model)
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/A{0}.npy'.format(i), A_list)
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/B{0}.npy'.format(i), B_list)
        model[model==0] = 1.
        flat = data/model
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat{0}.npy'.format(i), flat)
        plt.imshow(flat, vmin=0, vmax=1.5)
        plt.colorbar()
        plt.savefig('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat{0}.pdf'.format(i), dpi=190)
        plt.clf()
