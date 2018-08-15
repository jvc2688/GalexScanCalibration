import matplotlib
matplotlib.use('Agg')
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


if __name__ == '__main__':

#process
    if True:
        parser = argparse.ArgumentParser(description='galex scan count map')
        parser.add_argument('scan_name', nargs=1, help="scan name")
        parser.add_argument('np', nargs=1, type=int, help="number of process")
        parser.add_argument('pos', nargs=1, type=float, help="center of the image")
        parser.add_argument('out_path', nargs=1, help="path to the output directory")

        args = parser.parse_args()

        name_file = args.scan_name[0]

        num_t = args.np[0]
        skypos = [args.pos[0], 0.]

        out_path = args.out_path[0]


        print("scan name: {0}".format(name_file))
        print("skypos: {0}".format(skypos))
        print("output: {0}".format(out_path))
        
        pixsz = 0.0005555555556 #0.000416666666666667 #0.00166667
        skyrange = [10, 20]
        tranges = []
        trange = [0, 0]
        imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz).astype(int)
        wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
        
        count = np.zeros(imsz)

        for i in range(0, num_t):
            count += np.load('{0}/{1}_gal_sec_exp_tmp{2}.npy'.format(out_path, name_file, i)) #return_dict[i]
            #count += return_dict[i]

        tranges.append(trange)
        hdu = pyfits.PrimaryHDU(count)
        hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('{0}/count_map_{1}_exp.fits'.format(out_path, name_file), clobber=False)


        for i in range(num_t):
            os.remove('{0}/{1}_gal_sec_exp_tmp{2}.npy'.format(out_path, name_file, i))
            os.remove('/scratch/dw1519/galex/fits/scan_map/{1}_gal_sec_exp_tmp{2}.dat'.format(out_path, name_file, i))



