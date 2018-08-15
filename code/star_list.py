from __future__ import print_function
import catalog_fits as clog
from astropy.io import fits as pyfits
import numpy as np
import cPickle as pickle
import astropy.coordinates as pycoo
import os

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def get_star_pos(filename):
	data = np.loadtxt(filename, skiprows=4)
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

def get_star_pos_new(filename, limit):
	data = pyfits.open(filename)[1].data
	nuv = data['nuv_cntpsec']
	mask = (nuv<limit) & (nuv>20)
	data = data[mask]
	length = data.shape[0]
	stars = np.zeros((length,2))
	stars[:,0] = data['ra']
	stars[:,1] = data['dec']
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

def get_catalog_new(skypos, angle, stars_car, stars):
	center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
	rad = np.cos(angle*np.pi/180.)
	#print(center, rad)

	sep = np.dot(stars_car, center)
	coo = stars[sep>=rad, :]

	coo_list = []
	for i in range(0,coo.shape[0]):
		tmp_coo = (coo[i,0], coo[i,1])
		coo_list.append(tmp_coo)

	return coo_list


if __name__ == '__main__':
	if False:
		scan_list = ['AIS_GAL_SCAN_02183_0002','AIS_GAL_SCAN_02183_0003','AIS_GAL_SCAN_02192_0002','AIS_GAL_SCAN_02201_0001','AIS_GAL_SCAN_02201_0002','AIS_GAL_SCAN_02201_0003','AIS_GAL_SCAN_02210_0002']#['AIS_GAL_SCAN_00149_0001', 'AIS_GAL_SCAN_00761_0001', 'AIS_GAL_SCAN_01049_0001']#['05','14','23','32','41','50','59','68']
		name_file = 'name2147-2174'
		with open('../name_new/%s'%name_file) as f:
			scan_list = f.read().splitlines()
		print(scan_list)
		for scan_num in scan_list:

			asp_solution = np.load('../data/photon_list/%s_asp_cal_inter_half_sec.npy'%scan_num)

			stars_car, star_pos = get_star_pos('../data/sex_galex_matches_2147-2174_nofix.txt')

			star_set = set([])
			print(asp_solution.shape)
			for i in range(0, asp_solution.shape[0], 100):
				center = asp_solution[i,1:3]
				#stars = set(clog.get_catalog_t(center, 0.69))
				#stars = set(clog.get_catalog_t(center, 0.4))
				stars = set(get_catalog(center, 0.69, stars_car, star_pos))
				star_set.update(stars)
				print(i, len(star_set), center, len(stars))
			print(len(star_set))
			save_obj(star_set, '%s_starset_match'%scan_num)
			#save_obj(star_set, '%s_starset_extra_all_star'%scan_num)

	if True:
		finished = ['AIS_GAL_SCAN_02939_0002', 'AIS_GAL_SCAN_03209_0001', 'AIS_GAL_SCAN_03209_0002']

		with open('../check/finished', 'r') as f:
			finished = f.read().splitlines()

		problem = []
		with open('../check/problem', 'r') as f:
			problem = f.read().splitlines()

		scan_list = list(set(finished) - set(problem))
		#scan_list = ['AIS_GAL_SCAN_02489_0003', 'AIS_GAL_SCAN_02498_0003', 'AIS_GAL_SCAN_02588_0001']

		suffix = '-cal-sec'
		print(len(scan_list))

		for scan_num in scan_list:
			print(scan_num)
			output = '../data/{0}{1}-star/output.csv'.format(scan_num, suffix)
			dir = os.path.dirname(output)
			if not os.path.isfile(output):
				pass
			else:
				exit()
			if not os.path.exists(dir):
				os.makedirs(dir)
			else:
				print('exists')

			try:
				asp_solution = np.load('../data/photon_list/{0}{1}_asp_new.npy'.format(scan_num, suffix))
			except IOError:
				co_data = pyfits.open('../AIS_GAL_SCAN/asprta/{0}{1}-asprta.fits'.format(scan_num, suffix))[1].data
				T = co_data['T']
				ra = co_data['ra']
				dec = co_data['dec']
				roll = co_data['roll']
				t_new = np.arange((T.shape[0]-1)*200)*0.005+T[0]
				ra_new = np.interp(t_new, T, ra)
				dec_new = np.interp(t_new, T, dec)
				roll_new = np.interp(t_new, T, roll)
				asp_solution = np.array([t_new, ra_new, dec_new, roll_new]).T
				np.save('../data/photon_list/{0}{1}_asp_new.npy'.format(scan_num, suffix), asp_solution)

			#asp_solution = np.load('../data/photon_list/{0}{1}_asp_new.npy'.format(scan_num, suffix))#np.load('../data/photon_list/%s_asp_cal_inter_half_sec.npy'%scan_num)

			scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%scan_num)[1].data
			scst_time = scst_tmp['pktime']
			ix_tmp = np.digitize(asp_solution[:,0], scst_time) 
			ix_mask = ix_tmp<scst_time.shape[0]
			ix_tmp = ix_tmp[ix_mask]
			asp_solution = asp_solution[ix_mask]
			hv = scst_tmp['hvnom_nuv']
			mask = hv[ix_tmp]>0
			ix_tmp = ix_tmp[mask]
			asp_solution = asp_solution[mask]

			stars_car, star_pos = get_star_pos_new('../data/bstar.fits', 2000)

			star_set = set([])
			print(asp_solution.shape)
			for i in range(0, asp_solution.shape[0], 100):
				center = asp_solution[i,1:3]
				#stars = set(clog.get_catalog_t(center, 0.69))
				#stars = set(clog.get_catalog_t(center, 0.4))
				stars = set(get_catalog_new(center, 0.69, stars_car, star_pos))
				star_set.update(stars)
				print(i, len(star_set), center, len(stars))
			print(len(star_set))
			save_obj(star_set, '../data/{0}{1}-star/{0}_starset_new'.format(scan_num, suffix))
			#save_obj(star_set, '%s_starset_extra_all_star'%scan_num)

