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
from multiprocessing import Process
from multiprocessing import Manager
import re
import os

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

def update(pid, stars, exp, size, flat, return_dict):
    data = np.zeros((size,size))
    model = np.zeros((size, size))
    A_list = []
    for star in stars:
        print star
        star_mask = exp[:,0]==star
        p = exp[star_mask]
        x = p[:,1].astype(int)
        y = p[:,2].astype(int)
        A = np.sum(p[:,-2])/np.dot(p[:,-3],flat[x,y])
        A_list.append(A)
        model[x,y] += A*p[:,-3]
        data[x,y] += p[:,-2]    
    return_dict[pid] = (A_list, model, data)


if __name__ == '__main__':
    size=800
    num_t = 10
    outdir = sys.argv[1]#'super_iterate5'
    low = float(sys.argv[2])#0.015
    high = float(sys.argv[3])#0.02
    indir = 'super_exp'
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
    #flat[flat>0] = 1.
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
    #scans = ['0158', '0212', '0302', '0500', '0689', '0743', '1247', '1508', '1778', '2291', '2453', '2759', '2930' ,'3155', '3227', '3497']
    #scans = ['0158', '0212', '0302', '0500', '0689', '0743', '1247', '1508', '1778', '2291', '2453', '2759', '3227', '3497']
    
    scans = ['0005', '0023', '0050', '0158', '0212','0302', '0347', '0356', '0446', '0464',
             '0482', '0500', '0689', '0743','0806','0941', '0995','1247','1508', '1724',
             '1778', '1805', '2291', '2453', '2633', '2759', '2885', '3047', '3155', '3227', 
             '3326','3497']
    flist = glob.glob('/scratch/dw1519/galex/data/star_photon/super_exp/*_exp.npy')
    scans = [re.split('_', os.path.basename(f))[0] for f in flist]
    date = '08-21-2017'#'08-21-2017'
    #flat = np.load('/scratch/dw1519/galex/data/star_photon/iterate/flat3.npy')
    data = np.zeros((size,size))
    N=10
    for i in range(N):
        data = np.zeros((size,size))
        model = np.zeros((size,size))
        A_list = []
        for scan in scans:
            with open('../name_scan/%s'%scan) as f:
                name_list = f.read().splitlines()
            print name_list
            '''
            photon_list = []
            for name in name_list:
                photon_list.append(np.load('/scratch/dw1519/galex/data/star_photon/new/'+name+'_photon.npy'))
            photon_list = np.concatenate(photon_list, axis=0)
            '''
            #if i==0:
            #    data += np.load('/scratch/dw1519/galex/data/star_photon/full/'+scan+'_data.npy')
            exp = np.load('/scratch/dw1519/galex/data/star_photon/'+indir+'/'+scan+'_exp.npy')
            if len(exp.shape)<2:
                exp = np.concatenate(exp,axis=0)
            print exp.shape
            star_list = np.load('/scratch/dw1519/galex/data/star_photon/super/'+name_list[0]+'_star.npy')
            count_range = np.load('/scratch/dw1519/galex/data/star_photon/super/'+scan+'_mask.npy')
            count_rate = np.load('/scratch/dw1519/galex/data/star_photon/super/'+scan+'_ima.npy')

            print star_list.shape
            stars, ct = np.unique(exp[:,0], return_counts=True)#np.unique(photon_list[:,0], return_counts=True)
            nuv = star_list[stars.astype(int),-1]
            gb = star_list[stars.astype(int),-2]
            count_mask = np.ones(stars.shape[0], dtype=bool)
            count_list = np.zeros(stars.shape[0])-1.

            for j in range(count_rate.shape[0]):
                cr = count_rate[j]
                crange = count_range[j]
                gb_mask = (gb>=crange[0]) & (gb<=crange[1])
                count_list[gb_mask] = cr
            '''
            if count_range.shape[0]>0:
                for crange in count_range:
                    tmp_mask = ~((gb>=crange[0]) & (gb<=crange[1]))
                    count_mask &= tmp_mask
            '''
            count_mask = (count_list>=low) & (count_list<high)
            star_mask = (ct>10) &(nuv<17.5) &(nuv>15.5) & count_mask#(ct>10) & (nuv<18) &(nuv>16.)
            print np.sum(star_mask)
            '''
            gl = star_list[stars.astype(int),0]
            gb = star_list[stars.astype(int),1]
            bmax = np.max(gb)
            bmin = np.min(gb)
            lmax = np.max(gl)
            lmin = np.min(gl)
            if lmax-lmin > 2:
                gl[gl<350] = 360.-gl[gl<350]
                lmax = np.max(gl)
                lmin = np.min(gl)
            gl_mask = (gl>(lmin+0.3)) & (gl<(lmax-0.3))
            gb_mask = (gb>(bmin+0.5)) & (gb<(bmax-0.5))
            '''

            #print np.sum(star_mask)
            #stars = stars[star_mask&gl_mask&gb_mask]

            stars = stars[star_mask]
            
            stars_list = np.array_split(stars, num_t)
            manager = Manager()
            return_dict = manager.dict()
            p_list = []

            for pid in range(0, num_t):
                p = Process(target=update, args=(pid, stars_list[pid], exp, size, flat, return_dict))
                p.start()
                p_list.append(p)

            for p in p_list:
                p.join()

            for pid in range(0, num_t):
                result = return_dict[pid]
                A_list += result[0]
                model += result[1]
                data += result[2]      
        if i==0:
            np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/data.npy', data)
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/model{0}.npy'.format(i), model)
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/A{0}.npy'.format(i), A_list)
        model[model==0] = 1.
        flat = data/model
        np.save('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat{0}.npy'.format(i), flat)
        plt.imshow(flat, vmin=0, vmax=1.5)
        plt.colorbar()
        plt.savefig('/scratch/dw1519/galex/data/star_photon/'+outdir+'/flat{0}.pdf'.format(i), dpi=190)
        plt.clf()
