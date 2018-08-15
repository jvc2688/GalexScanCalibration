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
from astropy import units as u
import re

def determinant(v,w):
   return v[0]*w[1]-v[1]*w[0]

def ang(A, B):
    inner = np.arccos(np.dot(A, B))/np.pi*180.
    det = determinant(A,B)
    if det<0: #this is a property of the det. If the det < 0 then B is clockwise of A
        return inner
    else: # if the det > 0 then A is immediately clockwise of B
        return 360-inner

if __name__ == '__main__':

    asprta_path = '/scratch/dw1519/galex/AIS_GAL_SCAN/asprta'
    name = 'AIS_GAL_SCAN_01796_0003'
    asprta_suffix = '-cal-sec-dis-cut-correct-nuv-new'
    scst = '/scratch/dw1519/galex/AIS_GAL_SCAN/scst'
    scans = ['AIS_GAL_SCAN_00356_0001','AIS_GAL_SCAN_00392_0001','AIS_GAL_SCAN_00428_0001','AIS_GAL_SCAN_00437_0001','AIS_GAL_SCAN_00446_0001', 'AIS_GAL_SCAN_00455_0001', 'AIS_GAL_SCAN_00464_0001']
    scans = ['AIS_GAL_SCAN_00005_0001','AIS_GAL_SCAN_00014_0001','AIS_GAL_SCAN_00023_0001', 'AIS_GAL_SCAN_00032_0001', 'AIS_GAL_SCAN_00041_0001', 'AIS_GAL_SCAN_00050_0001',
            'AIS_GAL_SCAN_00059_0001', 'AIS_GAL_SCAN_00068_0001', 'AIS_GAL_SCAN_00086_0001', 'AIS_GAL_SCAN_00086_0002', 'AIS_GAL_SCAN_00095_0001',
            'AIS_GAL_SCAN_00095_0002', 'AIS_GAL_SCAN_00104_0001', 'AIS_GAL_SCAN_00104_0002', 'AIS_GAL_SCAN_03596_0001']
    scans = ['AIS_GAL_SCAN_00284_0001','AIS_GAL_SCAN_00293_0001','AIS_GAL_SCAN_00302_0001','AIS_GAL_SCAN_00311_0001','AIS_GAL_SCAN_00320_0001','AIS_GAL_SCAN_00329_0001','AIS_GAL_SCAN_00338_0001','AIS_GAL_SCAN_00347_0001', 'AIS_GAL_SCAN_00356_0001']
    scans = ['AIS_GAL_SCAN_00086_0001', 'AIS_GAL_SCAN_00095_0001', 'AIS_GAL_SCAN_00104_0001', 'AIS_GAL_SCAN_00113_0001',
            'AIS_GAL_SCAN_00122_0001', 'AIS_GAL_SCAN_00140_0001', 'AIS_GAL_SCAN_00149_0001', 'AIS_GAL_SCAN_00158_0001',
            'AIS_GAL_SCAN_00167_0001', 'AIS_GAL_SCAN_00176_0001', 'AIS_GAL_SCAN_00185_0001', 'AIS_GAL_SCAN_00194_0001']
    scans = ['AIS_GAL_SCAN_01724_0001', 'AIS_GAL_SCAN_01733_0002', 'AIS_GAL_SCAN_01742_0001', 'AIS_GAL_SCAN_01751_0001',
            'AIS_GAL_SCAN_01760_0001', 'AIS_GAL_SCAN_01769_0001', 'AIS_GAL_SCAN_01778_0001', 'AIS_GAL_SCAN_01787_0003',
            'AIS_GAL_SCAN_01796_0003', 'AIS_GAL_SCAN_01805_0002']

    scans = ['AIS_GAL_SCAN_02345_0002', 'AIS_GAL_SCAN_02354_0003', 'AIS_GAL_SCAN_02363_0001', 'AIS_GAL_SCAN_02372_0002',
            'AIS_GAL_SCAN_02381_0002', 'AIS_GAL_SCAN_02390_0002', 'AIS_GAL_SCAN_02399_0003', 'AIS_GAL_SCAN_02408_0002',
            'AIS_GAL_SCAN_02417_0002', 'AIS_GAL_SCAN_02426_0002', 'AIS_GAL_SCAN_02435_0002', 'AIS_GAL_SCAN_02444_0001']

    mean_gl = []
    mean_a = []
    for name in scans:
        co_data = pyfits.open('{0}/{1}{2}-asprta.fits'.format(asprta_path, name, asprta_suffix))[1].data
        T = co_data['T']
        ra = co_data['ra']
        dec = co_data['dec']
        roll = co_data['roll']
        asp_solution_tmp = np.array([T, ra, dec, roll]).T

        scst_tmp = pyfits.open('{0}/{1}-scst.fits'.format(scst, name))[1].data
        scst_time = scst_tmp['pktime']
        ix_tmp = np.digitize(asp_solution_tmp[:,0], scst_time)-1
        ix_mask = (ix_tmp>=0) & (ix_tmp<scst_time.shape[0])

        ix_tmp = ix_tmp[ix_mask]
        asp_solution_tmp = asp_solution_tmp[ix_mask]
        hv = scst_tmp['hvnom_nuv']
        mask = hv[ix_tmp]>0
        ix_tmp = ix_tmp[mask]
        asp_solution_tmp = asp_solution_tmp[mask]

        x = np.array([0,0])
        y = np.array([0,10000.])

        ra, dec = gn.gnomrev_simple(x[0], y[0], asp_solution_tmp[:,1], asp_solution_tmp[:,2], -asp_solution_tmp[:,3],1/36000.,0.)
        sky_data = SkyCoord(np.array([ra,dec]).T, unit='deg', frame=FK5, equinox='J2000.0')
        gal = sky_data.transform_to(Galactic)
        data0 = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)

        ra, dec = gn.gnomrev_simple(x[1], y[1], asp_solution_tmp[:,1], asp_solution_tmp[:,2], -asp_solution_tmp[:,3],1/36000.,0.)
        sky_data = SkyCoord(np.array([ra,dec]).T, unit='deg', frame=FK5, equinox='J2000.0')
        gal = sky_data.transform_to(Galactic)
        data1 = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
        up = np.array([0,1.])
        vector = ((data1-data0).T/np.linalg.norm((data1-data0),axis=1)).T
        angle = [ang(vector[i],up) for i in range(vector.shape[0])]
        mean_gl.append(np.mean(data0[:,0]))
        mean_a.append(np.mean(angle))

        plt.plot(angle, label=re.split('_',name)[3]+re.split('_',name)[4])
    plt.legend()
    plt.show()

    np.save('/scratch/dw1519/galex/fits/scan_map/11-29-2017/234-244.npy', np.array([mean_gl, mean_a]).T)
 

