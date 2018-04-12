import numpy as np
from astropy.io import fits as pyfits
import imagetools
import pos_range
import scipy
import math
import threading
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
import sys

if __name__ == '__main__':
	name = 'AIS_GAL_SCAN_02156_0002'
	name = sys.argv[1]
	print name
	asp = np.load('../data/photon_list/%s_asp_cal_inter_half_sec.npy'%name)[:,1:3]
	asp_o = np.load('../data/photon_list/%s_asp_cal_inter_half_sec_old.npy'%name)[:,1:3]
	gal = SkyCoord(asp, unit='deg', frame=FK5, equinox='J2000.0').transform_to(Galactic)
	gal_o = SkyCoord(asp_o, unit='deg', frame=FK5, equinox='J2000.0').transform_to(Galactic)
	asp_g = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
	asp_g_o = np.concatenate((np.array([gal_o.l.deg]).T, np.array([gal_o.b.deg]).T), axis=1)
	print np.mean(asp_g, axis=0)
	print np.mean(asp_g[:,0]-asp_g_o[:,0])*3600
	print np.mean(asp_g[:,1]-asp_g_o[:,1])*3600