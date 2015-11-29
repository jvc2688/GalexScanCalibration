from __future__ import print_function
import catalog_fits as clog
from astropy.io import fits as pyfits
import numpy as np
import cPickle as pickle
import astropy.coordinates as pycoo

def save_obj(obj, name ):
    with open('../data/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('../data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def get_star_pos(filename):
	data = np.loadtxt(filename, skiprows=1)
	print(data)
	length = data.shape[0]
	stars = np.zeros((length,5))
	stars[:,0] = data[:, 10]
	stars[:,1] = data[:, 11]
	stars[:,2] = data[:, 0]
	stars[:,3] = data[:, 3]
	stars[:,4] = data[:, 7]
	stars_rad = stars*np.pi/180.
	#convert (ra, dec) to (x, y, z)
	X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
	stars_car = np.array([X,Y,Z], dtype='float64').T
	return stars_car, stars

def get_catalog(skypos, angle, stars_car, stars):
	center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
	rad = np.cos(angle*np.pi/180.)
	#print(center, rad)

	sep = np.dot(stars_car, center)
	coo = stars[sep>=rad, :]

	coo_list = []
	for i in range(0,coo.shape[0]):
		tmp_coo = (coo[i,0], coo[i,1], coo[i,2], coo[i,3], coo[i,4])
		coo_list.append(tmp_coo)

	return coo_list


if __name__ == '__main__':
	scan_list = ['AIS_GAL_SCAN_00032_0001']#['AIS_GAL_SCAN_00149_0001', 'AIS_GAL_SCAN_00761_0001', 'AIS_GAL_SCAN_01049_0001']#['05','14','23','32','41','50','59','68']
	for scan_num in scan_list:

		asp_solution = np.load('../data/photon_list/%s_asp_cal.npy'%scan_num)

		stars_car, star_pos = get_star_pos('../data/sextractor_gl0-7.txt')

		star_set = set([])
		print(asp_solution.shape)
		for i in range(0, asp_solution.shape[0], 100):
			center = asp_solution[i,1:3]
			#stars = set(clog.get_catalog_t(center, 0.69))
			#stars = set(clog.get_catalog_t(center, 0.4))
			stars = set(get_catalog(center, 0.4, stars_car, star_pos))
			star_set.update(stars)
			print(i, len(star_set))
		print(len(star_set))
		save_obj(star_set, '%s_starset_extra_full'%scan_num)

