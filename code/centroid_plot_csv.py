#import matplotlib
#matplotlib.use('Agg')
import c3
from astropy.io import fits as pyfits
import numpy as np
from astropy import wcs as pywcs
import matplotlib.pyplot as plt
import sys
import aplpy
import os
from sklearn.neighbors import KernelDensity
import csv
import math
import asp_cal
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude

def _find_centroid(filename):
	try:
		hdulist = pyfits.open(filename)
	except IOError as e:
		print "I/O error({0}): {1}".format(e.errno, e.strerror)
		hdulist = None
	except:
		print "Unexpected error:", sys.exc_info()[0]
 		raise
	if hdulist is not None:
		w = pywcs.WCS(hdulist[0].header, hdulist)
		data = hdulist[0].data
		data = data.byteswap(True).newbyteorder()
		cy, cx = c3.find_centroid(data)
	else:
		return None
	return cx,cy

def find_centroid(filename):
	try:
		hdulist = pyfits.open(filename)
	except IOError as e:
		print "I/O error({0}): {1}: {2}".format(e.errno, e.strerror, filename)
		hdulist = None
	except:
		print "Unexpected error:", sys.exc_info()[0]
		raise
	if hdulist is not None:
		hdulist = pyfits.open(filename)
		w = pywcs.WCS(hdulist[0].header, hdulist)
		data = hdulist[0].data
		data = data.byteswap(True).newbyteorder()
		cy, cx = c3.find_centroid(data)
		centroid = w.wcs_pix2world(w.sip_foc2pix([[cx, cy]],1),1)[0]
		if centroid[0]>1:
			centroid[0] = centroid[0]-360.
	else:
		centroid = [0,0]
	return centroid

def get_centers(initial, final):
	centers = []
	for i in range(initial, final+1):
		filename = '../fits/co/right/co_map%d_%d_zoom_nocal.fits'%(i,i+1)
		centroid = find_centroid(filename)
		centers.append(centroid)
	centers = np.array(centers)
	return centers

def get_centers_half(initial, final):
	centers = []
	for i in range(initial, final+1):
		for j in range(10):
			filename = '../fits/test/co_map%d_%d_10.fits'%(i,j)
			centroid = find_centroid(filename)
			centers.append(centroid)
	centers = np.array(centers)
	return centers

def corr_plot(centroid, filename, title):

	print filename
	#centroid = find_centroid(filename)
	#centroid = np.load('../data/offsets300-899_half_r.npy')[7]
	print centroid
	fig = aplpy.FITSFigure(filename)
	fig.add_label(centroid[0], centroid[1], 'X', color='red')
	fig.show_grayscale(invert=True)
	fig.tick_labels.set_xformat('d.ddd')
	fig.tick_labels.set_yformat('d.ddd')
	fig.recenter(0., 0., radius=0.01)
	fig.add_grid()
	fig.set_title(title)
	basename = os.path.basename(filename)
	preffix, ext = os.path.splitext(basename)
	fig.save('../plots/corr_10/%s.png'%preffix)

def load_data(filename):
	data = []
	with open(filename, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			if math.fabs(float(row[0]))>0.008 or math.fabs(float(row[1]))>0.008:
				pass
			else:
				data.append(row)
		csvfile.close()
	return data

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

def generate_zero_offsets(name):
	print name
	hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
	initial = 1
	final = hdulist[1].data['T'].shape[0]-1

	centers = []
	for i in range(initial, final+1):
		centers.append(np.load('../data/%s/cata/centroids_photon%d.npy'%(name, i)))
	centroids = np.concatenate(centers, axis=0)

	print centroids.shape

	np.save('../data/%s/cata/offsets%d_10_new_photon.npy'%(name, initial), centroids)

	asp_cal.interpolate_offsets(name)

	output = "../plots/%s/cata/output.csv"%(name)
	dir = os.path.dirname(output)
	if not os.path.exists(dir):
		os.makedirs(dir)


def generate_first_offsets(name):
	print name
	hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
	initial = 1
	final = hdulist[1].data['T'].shape[0]-1

	centers = []
	for i in range(initial, final+1):
		centers.append(np.load('../data/%s/cata/centroids%d_half.npy'%(name, i)))
	centroids = np.concatenate(centers, axis=0)

	print centroids.shape

	out_mask = np.zeros(centroids.shape[0])

	z, mz = moving_stat(centroids[:,0], out_mask, half_win=100)
	outliers = np.zeros(centroids.shape[0])
	outliers[mz>3.5] = 1
	outliers[out_mask>0] = 1
	outliers = outliers>0
	index = np.arange(centroids.shape[0])
	centroids[outliers, 0] = np.interp(index[outliers], index[~outliers], centroids[~outliers,0])

	z, mz = moving_stat(centroids[:,1], out_mask, half_win=100)
	outliers = np.zeros(centroids.shape[0])
	outliers[mz>3.5] = 1
	outliers[out_mask>0] = 1
	outliers = outliers>0
	index = np.arange(centroids.shape[0])
	centroids[outliers, 1] = np.interp(index[outliers], index[~outliers], centroids[~outliers,1])

	plt.plot(centroids[:,0], '.b')
	plt.show()

	np.save('../data/%s/cata/offsets%d_10_new_half.npy'%(name, initial), centroids)

	output = "../plots/%s/cata/output.csv"%(name)
	dir = os.path.dirname(output)
	if not os.path.exists(dir):
		os.makedirs(dir)

	asp_cal.interpolate_offsets(name, 0.5, centroids)


def generate_sec_offsets(name):
	print name
	hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
	initial = 1
	final = hdulist[1].data['T'].shape[0]-1

	centers = []
	for i in range(initial, final+1):
		centers.append(np.load('../data/%s/cata/centroids_sec%d.npy'%(name, i)))
	centroids = np.concatenate(centers, axis=0)

	print centroids.shape
	plt.plot(centroids[:,0], '.b')
	plt.show()

	np.save('../data/%s/cata/offsets%d_10_new_sec.npy'%(name, initial), centroids)

	asp_cal.secondary_cal(name)

if __name__ == '__main__':
	if False:
		name = 'AIS_GAL_SCAN_00104_0001'
		generate_first_offsets(name)

	if True:
		name = 'AIS_GAL_SCAN_00104_0001'
		generate_sec_offsets(name)

	if False:
		'''
		scan_num = sys.argv[1]
		initial = int(sys.argv[2])#300
		final = int(sys.argv[3])#1343
		'''
		name_file = sys.argv[1]
		print name_file
		with open('../name/%s'%name_file) as f:
			name_list = f.read().splitlines()

		for name in name_list:
			print name
			hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
			initial = 1
			final = hdulist[1].data['T'].shape[0]-1

			centers = []
			for i in range(initial, final+1):
				centers.append(np.load('../data/%s/cata/centroids%d.npy'%(name, i)))
			centroids = np.concatenate(centers, axis=0)
			#centroids = -centroids
			print centroids.shape

			trange = np.arange((final-initial+1)*10)/10.+initial-0.5
			np.save('../data/%s/cata/offsets%d_10_new.npy'%(name, initial), centroids)

			asp_cal.run_cal(name)

			output = "../plots/%s/cata/output.csv"%(name)
			dir = os.path.dirname(output)
			if not os.path.exists(dir):
				os.makedirs(dir)
			
			plt.plot(trange, centroids[:,0],'.b', markersize=2)
			plt.xlabel('time/s')
			plt.xlim(initial-2,final+2)
			#plt.ylim(-0.005, 0.005)
			plt.ylabel('RA/degree')
			plt.savefig('../plots/%s/cata/ra%d_%d_new_r.png'%(name,initial,final),dpi=190)
			plt.clf()
			plt.plot(trange, centroids[:,1],'.b', markersize=2)
			plt.xlabel('time/s')
			plt.xlim(initial-2,final+2)
			#plt.ylim(-0.005, 0.005)
			plt.ylabel('DEC/degree')
			plt.savefig('../plots/%s/cata/dec%d_%d_new_r.png'%(name,initial,final),dpi=190)
			plt.clf()

#indentify outliers
	if False:
		'''
		scan_num = sys.argv[1]
		initial = int(sys.argv[2])#300
		final = int(sys.argv[3])#1343
		'''
		name_file = sys.argv[1]
		print name_file
		with open('../name_new/%s'%name_file) as f:
			name_list = f.read().splitlines()

		for name in name_list:
			print name
			hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
			initial = 1
			final = hdulist[1].data['T'].shape[0]-1

			centers = []
			masks = []

			for i in range(initial, final+1):
				centers.append(np.load('../data/%s/cata/centroids_sec%d.npy'%(name, i)))
				#masks.append(np.load('../data/%s/cata/out_mask%d.npy'%(name, i)))
			centroids = np.concatenate(centers, axis=0)

			#centroids = np.load('../data/%s/cata/offsets%d_10_new_inter_half_fine.npy'%(name, initial))

			#centroids_fine = np.load('../data/%s/cata/offsets%d_10_new.npy'%(name, initial))
			#out_mask = np.concatenate(masks, axis=0)
			out_mask = np.zeros(centroids.shape[0])
			'''
			z, mz = moving_stat(centroids[:,0], out_mask, half_win=250)
			outliers = np.zeros(centroids.shape[0])
			outliers[mz>3.5] = 2
			outliers[out_mask>0] = 2
			outliers = outliers>1
			index = np.arange(centroids.shape[0])
			centroids[outliers, 0] = np.interp(index[outliers], index[~outliers], centroids[~outliers,0])
			'''

			'''
			z, mz = moving_stat(centroids[:,0], out_mask, half_win=150)
			outliers_l = np.zeros(centroids.shape[0])
			outliers_l[mz>3.5] = 1
			outliers_l[out_mask>0] = 1

			outliers = outliers+outliers_l
			'''

			'''
			z, mz = moving_stat(centroids[:,1], out_mask, half_win=250)
			outliers_l = np.zeros(centroids.shape[0])
			outliers_l[mz>3.5] = 2
			outliers_l[out_mask>0] = 2
			outliers_l = outliers_l>1
			centroids[outliers_l, 1] = np.interp(index[outliers_l], index[~outliers_l], centroids[~outliers_l,1])
			'''

			#centroids = -centroids
			print centroids.shape

			trange_fine = np.arange((final-initial+1)*10)/10.+initial-0.5
			np.save('../data/%s/cata/offsets%d_10_new_sec.npy'%(name, initial), centroids)

			asp_cal.secondary_cal(name)

			output = "../plots/%s/cata/output.csv"%(name)
			dir = os.path.dirname(output)
			if not os.path.exists(dir):
				os.makedirs(dir)
			
			trange = np.arange(centroids.shape[0], dtype=float) + initial-0.5
			f, axes = plt.subplots(2, 1, squeeze=False)
			axes[0,0].plot(trange, centroids[:,0],'.b', markersize=2)
			#axes[0,0].plot(trange_fine, centroids_fine[:,0],'.r', markersize=2)

			#axes[0,0].plot(trange[outliers>1], centroids[:,0][outliers>1],'.r', markersize=2)
			#axes[0,0].plot(trange[outliers_z>0], centroids[:,0][outliers_z>0],'.g', markersize=4)

			plt.xlabel('time/s')
			plt.xlim(initial-2,final+2)
			#plt.ylim(-0.005, 0.005)
			plt.ylabel('RA/degree')
			#plt.savefig('../plots/%s/cata/ra%d_%d_new_o.png'%(name,initial,final),dpi=190)
			axes[1,0].plot(trange, centroids[:,1],'.b', markersize=2)
			#axes[1,0].plot(trange_fine, centroids_fine[:,1],'.r', markersize=2)

			#axes[1,0].plot(trange[outliers_l>1], centroids[:,1][outliers_l>1],'.r', markersize=2)
			plt.xlabel('time/s')
			plt.xlim(initial-2,final+2)
			#plt.ylim(-0.005, 0.005)
			plt.ylabel('DEC/degree')
			#plt.savefig('../plots/%s/cata/dec%d_%d_new_o.png'%(name,initial,final),dpi=190)
			plt.show()
		
	if False:
		initial = 800
		final = 800
		centers = np.load('../data/offsets300-1343_10_new.npy')
		for i in range(initial, final+1):
			for j in range(10):
				filename = '../fits/test/co_map%d_%d_10.fits'%(i,j)
				corr_plot(centers[(i-300)*10+j],filename, 'corr between %.1fs and %.1fs'%(i+0.1*j-0.5,i+0.1*(j+1)-0.5))

	if False:
		initial = 300
		final = 1343

		trange = np.arange((final-initial+1)*10)/10.+initial-0.5
		centers = np.load('../data/05_r/cata/offsets%d-%d_10_new.npy'%(initial,final))

		'''
		for i in range(1, (final-initial+1)*10):
			centers[i] = centers[i]+centers[i-1]
		'''
		#centers = -centers

		centers = np.diff(centers, axis=0)

		plt.clf()
		plt.plot(trange, centers[:,0],'.b', markersize=2)
		plt.xlabel('time/s')
		plt.xlim(initial-2, final+2)
		plt.ylim(-0.03,0.12)
		plt.ylabel('RA/degree')
		plt.savefig('../plots/photon_05_r/cata/ra%d_%d_cul.png'%(initial,final),dpi=190)
		plt.clf()
		plt.plot(trange, centers[:,1],'.b', markersize=2)
		plt.xlabel('time/s')
		plt.xlim(initial-2, final+2)
		plt.ylim(-0.02, 0.05)
		plt.ylabel('DEC/degree')
		plt.savefig('../plots/photon_05_r/cata/dec%d_%d_cul.png'%(initial,final),dpi=190)

	if False:

		scan_name = 'AIS_GAL_SCAN_00149_0001'
		asp = np.load('../data/photon_list/%s_asp_cal_inter_half.npy'%scan_name)
		asp_uni, uni_index=np.unique(asp[:,0], True)
		asp_uni = asp[uni_index]
		asp_o = np.load('../data/photon_list/%s_asp.npy'%scan_name)
		asp_uni_o, uni_index_o=np.unique(asp_o[:,0], True)
		asp_uni_o = asp_o[uni_index_o]

		#asp_diff = np.diff(asp_uni[:,1:3], axis=0)
		'''
		sky_data = SkyCoord(asp_uni[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
		gal = sky_data.transform_to(Galactic)
		asp_uni[:,1:3] = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1) 

		sky_data = SkyCoord(asp_uni_o[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
		gal = sky_data.transform_to(Galactic)
		asp_uni_o[:,1:3] = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1) 
		'''
		f, axes = plt.subplots(2, 1, squeeze=False)
		axes[0,0].plot(asp_uni[:,0]-asp_uni[0,0], asp_uni[:,1],'.r')
		axes[1,0].plot(asp_uni[:,0]-asp_uni[0,0], asp_uni[:,2],'.r')
		axes[0,0].plot(asp_uni_o[:,0]-asp_uni_o[0,0], asp_uni_o[:,1],'.k')
		axes[1,0].plot(asp_uni_o[:,0]-asp_uni_o[0,0], asp_uni_o[:,2],'.k')
		#plt.setp( axes[0,0].get_xticklabels(), visible=False)
		plt.show()
