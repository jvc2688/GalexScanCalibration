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

asp_list = ['05','14','23','32','41','50','59','68']

ra_min, ra_max, dec_min, dec_max = [],[],[],[]
for asp_num in asp_list:
	asp = np.load('../data/photon_list/NUVPhoton%s_full_new_asp_cal_cata_wa.npy'%asp_num)
	ra_min.append(np.min(asp[:,1]))
	ra_max.append(np.max(asp[:,1]))
	dec_min.append(np.min(asp[:,2]))
	dec_max.append(np.max(asp[:,2]))

print np.min(ra_min), np.max(ra_max), np.min(dec_min), np.max(dec_max)
print (np.min(ra_min)+np.max(ra_max))/2, (np.min(dec_min)+np.max(dec_max))/2
print np.max(ra_max)-np.min(ra_min), np.max(dec_max)-np.min(dec_min)