import co_rel_csv as crc
from sklearn.neighbors.kde import KernelDensity
import numpy as np
from sklearn.cross_validation import KFold
from scipy import stats
import sys
from astropy.io import fits as pyfits

if __name__ == '__main__':
	initial_sec = int(sys.argv[1])
	step = 0.1
	num_co = int(1/step)

	hdulist = pyfits.open('../data/AIS_GAL_SCAN_00005_0001-asprta.fits')
	co_data = hdulist[1].data
	intitial_asp = co_data[initial_sec]
	center = np.array([intitial_asp[1], intitial_asp[2]])  

	skypos = [0.0, 0.0]
	skyrange = [0.04, 0.04]

	initial_time = intitial_asp[0]-0.5

	tranges = []
	for sec in range(num_co+1):
		tranges.append([initial_time+step*sec, initial_time+step*(sec+1)])
	print tranges

	data = crc.get_data(tranges)
	print len(data)

  	centroids = []
  	for sec in range(num_co):
		coo1 = np.array(data[sec], dtype='float64')[:,-3:-1]
		coo2 = np.array(data[sec+1], dtype='float64')[:,-3:-1]

		coo1 = crc.angle_filter(coo1, center, 1.)
		coo2 = crc.angle_filter(coo2, center, 1.)

		co_rel = crc.get_corr_map_(coo2, coo1, skypos, skyrange, sec)
    	
		kernel = stats.gaussian_kde(co_rel)

		print kernel.factor, kernel.d, kernel.n

