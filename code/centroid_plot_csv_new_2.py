import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import c3
from astropy.io import fits as pyfits
import numpy as np
from astropy import wcs as pywcs
import sys
import aplpy
import os
from sklearn.neighbors import KernelDensity
import csv
import math
import asp_cal
import glob
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.io.fits import update
import re
from scipy.interpolate import splev, splrep


def moving_stat(data, out_mask, half_win=10):
    moving_mask = np.zeros(data.shape)
    moving_mask[0:2*half_win+1] = 1
    mean = np.zeros(data.shape)
    median = np.zeros(data.shape)
    abs_dev = np.zeros(data.shape)
    std = np.zeros(data.shape)
    z = np.zeros(data.shape)
    mz = np.zeros(data.shape)

    for i in range(half_win):
        if out_mask[i] == 1:
            std[i] = 1
            abs_dev[i] = 1
        else:
            tmp_out_mask = -out_mask[:half_win+i+1]+1
            #print i, data[:half_win+i+1][tmp_out_mask>0].shape
            mean[i] = np.mean(data[:half_win+i+1][tmp_out_mask>0], axis=0)
            std[i] = np.std(data[:half_win+i+1][tmp_out_mask>0], axis=0)

            median[i] = np.median(data[:half_win+i+1][tmp_out_mask>0], axis=0)
            abs_dev[i] = np.median(np.absolute(data[:half_win+i+1][tmp_out_mask>0]-median[i]), axis=0)

        if out_mask[-i-1] == 1:
            std[-i-1] = 1
            abs_dev[-i-1] =1
        else:
            tmp_out_mask = -out_mask[-half_win-i-1:]+1
            median[-i-1] = np.median(data[-half_win-i-1:][tmp_out_mask>0], axis=0)
            abs_dev[-i-1] = np.median(np.absolute(data[-half_win-i-1:][tmp_out_mask>0]-median[-i-1]), axis=0)

            mean[-i-1] = np.mean(data[-half_win-i-1:][tmp_out_mask>0], axis=0)
            std[-i-1] = np.std(data[-half_win-i-1:][tmp_out_mask>0], axis=0)
            #print -i-1, data[-half_win-i-1:][tmp_out_mask>0].shape

    for i in range(data.shape[0]-2*half_win):
        if out_mask[half_win+i] == 1:
            std[half_win+i] = 1
            abs_dev[half_win+i] =1
        moving_mask_tmp = np.roll(moving_mask, i)
        tmp_out_mask = -out_mask[moving_mask_tmp>0]+1
        mean[half_win+i] = np.mean(data[moving_mask_tmp>0][tmp_out_mask>0], axis=0)
        std[half_win+i] = np.std(data[moving_mask_tmp>0][tmp_out_mask>0], axis=0)
        median[half_win+i] = np.median(data[moving_mask_tmp>0][tmp_out_mask>0], axis=0)
        abs_dev[half_win+i] = np.median(np.absolute(data[moving_mask_tmp>0][tmp_out_mask>0]-median[half_win+i]), axis=0)

        #print half_win+i, data[moving_mask_tmp>0][tmp_out_mask>0].shape

    z = np.absolute((data - mean)/std)
    mz = np.absolute(0.6745*(data-median)/abs_dev)

    return z, mz


def generate_new_offsets_new(name, asprta, suffix, tmp_dir, num_p, hv_mask, nan_mask):
    print name

    try:
        centroids = np.load(tmp_dir+'/offsets_%s.npy'%(suffix))
        time = np.load(tmp_dir+'/time_%s.npy'%(suffix))
    except IOError:
        try:
            centroids = np.load(tmp_dir+'/offsets0_%s.npy'%(suffix))
            time = np.load(tmp_dir+'/time0_%s.npy'%(suffix))
        except IOError:
            print 'no file'
            return 0

    print centroids.shape
    print time.shape

    output = '../plots/2/%s-%s.pdf'%(name, suffix)
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
      os.makedirs(dir)

    print np.sum(nan_mask)
    if np.sum(nan_mask)>0:
        x = np.arange(centroids.shape[0])
        ct = centroids.copy()
        ct[nan_mask] = np.mean(ct[hv_mask][~nan_mask[hv_mask]], axis=0)
        plt.plot(x,centroids[:,0], '.k')
        plt.plot(x[nan_mask],ct[nan_mask,0], '.r')
        plt.savefig('../plots/2/%s-%s_ra.pdf'%(name, suffix), dpi=190)
        plt.clf()

        plt.plot(centroids[:,1], '.k')
        plt.plot(x[nan_mask], ct[nan_mask,1], '.r')
        plt.savefig('../plots/2/%s-%s_dec.pdf'%(name, suffix), dpi=190)
        plt.clf()
        centroids = ct
    else:
        plt.plot(centroids[:,0], '.k')
        plt.savefig('../plots/2/%s-%s_ra.pdf'%(name, suffix), dpi=190)
        plt.clf()

        plt.plot(centroids[:,1], '.k')
        plt.savefig('../plots/2/%s-%s_dec.pdf'%(name, suffix), dpi=190)
        plt.clf()

    #np.save(tmp_dir+'/offsets_%s.npy'%(suffix), centroids)
    #np.save(tmp_dir+'/time_%s.npy'%(suffix), time)
    
    hdulist = pyfits.open(asprta)
    co_data = hdulist[1].data
    T = co_data['T']
    ra = co_data['ra']
    dec = co_data['dec']
    roll = co_data['roll']
    ra_new = np.interp(time, T, ra) - centroids[:,0]
    dec_new = np.interp(time, T, dec) - centroids[:,1]
    roll_new = np.interp(time, T, roll) - centroids[:,2]

    other = np.zeros((time.shape[0], 8))
    array = np.concatenate([np.array([time, ra_new, dec_new, roll_new]).T, other], axis=1)
    data = np.core.records.fromarrays(array.transpose(), dtype=[('T', float), ('RA', float), ('DEC', float), ('ROLL', float),\
                 ('STATUS_FLAG', int), ('ROLL_RAD', float), ('X', float), ('Y', float), ('Z', float), ('XDOT', float), ('YDOT', float), ('ZDOT', float)])

    new_file = re.split('-asprta.fits', asprta)[0]+'-'+suffix+'-asprta.fits'
    os.system('cp {0} {1}'.format(asprta, new_file))
    update(new_file, data, 1)
    hdu = pyfits.open(new_file)
    print hdu[1].data['RA'].shape
    print hdu[1].data['DEC'].shape
    hdu.close()



    
if __name__ == '__main__':
    if True:
        name = sys.argv[1]
        suffix = sys.argv[2]
        generate_new_offsets_new(name, suffix)
