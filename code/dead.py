import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits
import imagetools
import gnomonic as gn
import pos_range
import scipy
import math
import threading
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from multiprocessing import Process
from multiprocessing import Manager
import sys
import os
import gc
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
import argparse


if __name__ == '__main__':
    sub_name = sys.argv[1]
    print sub_name
    scan_list = []
    with open(sub_name, 'r') as f:
        scan_list = f.read().splitlines()
    for name in scan_list:
        print name
        scst = '/scratch/dw1519/galex/AIS_GAL_SCAN/scst'
        scst_tmp = pyfits.open('{0}/{1}-scst.fits'.format(scst, name))[1].data
        scst_time = scst_tmp['pktime']

        hv = scst_tmp['hvnom_nuv']

        dead = scst_tmp['t_dead_nuv']
        tec = scst_tmp['NDCTEC']
        fec = scst_tmp['NDCFEC']
        plt.plot(fec, 1-dead, '.k', label='1-dead')
        plt.plot(fec, tec/fec, '.r', label='TEC/FEC')
        plt.legend()
        plt.savefig('/scratch/dw1519/galex/plots/dead/'+name+'.png', dpi=190)
        plt.clf()

        hold = 150000
        fec_mask = fec>hold
        width=1000
        while True:
            ratio_mask = (fec>(hold-width)) & (fec<(hold+width))
            if np.sum(ratio_mask)>10:
                break
            else:
                width+=1000
        ratio = np.mean((1.-dead[ratio_mask])/(tec[ratio_mask]/fec[ratio_mask]))
        print 'ratio:{0}'.format(ratio)
        ratio = np.median(1.-dead[ratio_mask])/np.median(tec[ratio_mask]/fec[ratio_mask])
        print 'ratio:{0}'.format(ratio)
        #dead[fec_mask] = 1.-ratio*tec[fec_mask]/fec[fec_mask]
        plt.plot(scst_time, dead, '.k')
        plt.plot(scst_time[fec_mask], 1.-ratio*tec[fec_mask]/fec[fec_mask], '.r')
        plt.savefig('/scratch/dw1519/galex/plots/dead/'+name+'_dead.png', dpi=190)
        plt.clf()


