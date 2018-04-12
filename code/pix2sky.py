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

if __name__ == '__main__':

  if True:
    filename = '../data/sex_total_203-257.txt'
    data = np.loadtxt(filename, skiprows=1)
    pix_coo = data[:,1:3]

    name_file = sys.argv[1]
    with open('../name_new/%s'%name_file) as f:
      name_list = f.read().splitlines()
    print name_list

    csv_lists = {}
    total_num_csv = 0
    gal_l = []
    for name in name_list:
      csvs = glob.glob("../data/%s/split/*.csv"%name)
      csv_lists[name] = csvs
      total_num_csv += len(csvs)

      asprta = np.load('../data/photon_list/%s_asp_cal.npy'%name)#np.load('../data/photon_list/%s_asp.npy'%name)
      sky_data = SkyCoord(asprta[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
      gal = sky_data.transform_to(Galactic)
      asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
      gal_l.append(np.mean(asprta[:,0]))

    skypos = [np.mean(gal_l), 0.]
    skyrange = [np.max(gal_l)-np.min(gal_l)+2., 24.]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix_g(skypos, skyrange, 0.00166667)
    wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.00166667)

    sky_coo = wcs.wcs_pix2world(wcs.sip_foc2pix(pix_coo,1),1)
    print sky_coo
