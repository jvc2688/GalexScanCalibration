import sys
import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc
from astropy.io import fits as pyfits
import catalog_fits
import astropy.coordinates as pycoo
import centroid_kde as ck
import cPickle as pickle
import catalog_fits
import os
from multiprocessing import Process
from multiprocessing import Manager
import spilt_csv as spi
import gnomonic as gn
import leastSquareSolver as lss
from scipy import ndimage

def split_seq(seq, size):
  newseq = []
  splitsize = 1.0/size*len(seq)
  for i in range(size):
    newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
  return newseq

def get_detector_pos(pos_list):
	pos_array = np.array(pos_list, dtype='float64')
	pos = ((pos_array/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800
	return pos

def get_res(pos_list):
	flat = pyfits.open('../data/cal/NUV_flat.fits')[0].data
	flat = np.swapaxes(flat,0,1)
	pos_array = np.array(pos_list, dtype='float64')
	pos = np.mean(pos_array, axis=0)
	pos = ((pos/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800
	return flat[pos[0], pos[1]]

def get_res_new(xi, eta, file_name='../data/cal/NUV_flat.fits'):
	flat = pyfits.open(file_name)[0].data
	#flat = np.swapaxes(flat,0,1)
	x = ((xi/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800
	y = ((eta/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800
	res = np.zeros(x.shape[0])
	for i in range(x.shape[0]):
		for l in range(-1,2):
			for m in range(-1,2):
				x0 = x[i]+l
				y0 = y[i]+m
				if x0>=800 or x0<0 or y0>=800 or y0<0:
					return None
				res[i] += flat[x0, y0]
	return res

def get_pos(pos_list):
	pos_array = np.array(pos_list, dtype='float64')
	pos = np.mean(pos_array, axis=0)
	pos = ((pos/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800
	return pos

def res_star(pid, scan_name, start, end, return_dict):
	cata = spi.load_obj('../data/%s_starset_new'%name)
	cata_a = np.array(list(cata))
	cata_len = len(cata)
	cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)
  
	for star_num in range(start, end+1): 
		order = 0
		star = 'star_%d-%d'%(star_num, order)
		
		if not os.path.exists('../data/%s/star/list/%s.csv'%(scan_name, star)):
			print 'skip:%s'%star
			continue
	  
		csv_file = '../data/%s/star/list/%s.csv'%(scan_name, star)
		data = []
		time_list = []
		with open(csv_file, 'rb') as file:
			reader = csv.reader(file)
			first = reader.next()
			last = float(first[0])
			time_list.append(last)
			data.append([])
			data[0].append(first[6:8])
			i = 0

			for row in reader:
			  time = float(row[0])
			  if time - last > 500:
			    data.append([])
			    i += 1
			    data[i].append(row[6:8])
			    last = time
			    time_list.append(last)
			  else:
			    data[i].append(row[6:8])
		
		output = "../plots/%s/star/res/%s.csv"%(scan_name, scan_name)
		dir = os.path.dirname(output)
		if not os.path.exists(dir):
			os.makedirs(dir)

		res_list = map(get_res, data)
		print map(get_pos, data)
		'''
		plt.plot(time_list, res_list)
		plt.xlabel('t/s')
		plt.ylabel('response')
		plt.savefig('../plots/%s/star/res/res%d-%d.png'%(scan_name, star_num, order), dpi=190)
		plt.clf()
		'''

def res_star_new(pid, scan_name, start, end, return_dict):
	cata = spi.load_obj('../data/%s_starset_new'%scan_name)
	cata_a = np.array(list(cata))
	cata_len = len(cata)
	cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)

	asp = np.load('../data/photon_list/AIS_GAL_SCAN_00014_0001_asp_cal.npy')
	asp_uni, uni_index=np.unique(asp[:,0], True)
	asp_uni = asp[uni_index]
	ix_tmp = (np.round(asp_uni[:,0]-asp_uni[0,0])+1).astype(int)
	dead_tmp = np.zeros(ix_tmp.shape[0])
	scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%scan_name)[1].data
	limit = scst_tmp['t_dead_nuv'].shape[0]-1
	ix_tmp[ix_tmp>limit] = limit
	dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]

	for star_num in range(start, end+1): 
		star_co = cata_a[star_num, 0:2]


		order = 1
		star = 'star_%d-%d'%(star_num, order)
		
		if not os.path.exists('../data/%s/star/list/%s.csv'%(scan_name, star)):
			print 'skip:%s'%star
			continue
	  
		csv_file = '../data/%s/star/list/%s.csv'%(scan_name, star)
		count = []
		time_list = []
		with open(csv_file, 'rb') as file:
			reader = csv.reader(file)
			first = float(reader.next()[0])
			final = 0
			last = 0
			reader = csv.reader(file)
			for row in reader:
				time = float(row[0])#time = int(float(row[0])/1000.)
				if time - last > 200:
					count.append(1)
					last = time
					time_list.append(last)
				else:
					count[-1] += 1
				final = time

		time_list = np.array(time_list)/1000.
		index_0 = (first-asp_uni[0,0]*1000)/5
		index_1 = (final-asp_uni[0,0]*1000)/5

		star_track = asp_uni[index_0:index_1+1,:]
		dead_t = dead_tmp[index_0:index_1+1]
		ra = np.ones(star_track.shape[0])*star_co[0]
		dec = np.ones(star_track.shape[0])*star_co[1]

		xi, eta = gn.gnomfwd_simple(ra, dec, star_track[:,1], star_track[:,2], -star_track[:,3], 1/36000., 0.0)

		output = "../plots/%s/star/res/%s.csv"%(scan_name, scan_name)
		dir = os.path.dirname(output)
		if not os.path.exists(dir):
			os.makedirs(dir)
		print dead_t
		res_list = get_res_new(xi, eta)
		res_list_d = res_list*(1-dead_t)

		'''
		plt.plot(star_track[:,0], res_list_d)
		plt.xlabel('t/s')
		plt.ylabel('response')
		plt.savefig('../plots/%s/star/res/res%d-%d_track_d.png'%(scan_name, star_num, order), dpi=190)
		plt.clf()

		plt.plot(star_track[:,0], res_list)
		plt.xlabel('t/s')
		plt.ylabel('response')
		plt.savefig('../plots/%s/star/res/res%d-%d_track.png'%(scan_name, star_num, order), dpi=190)
		plt.clf()
		'''

		f, axes = plt.subplots(2, 1)
		axes[0].plot(time_list, count)
		axes[0].set_ylabel('Num of Photons')
		ylim = axes[0].get_ylim()
		axes[0].set_ylim(ylim+np.absolute(ylim)*0.01)
		plt.setp( axes[0].get_xticklabels(), visible=False)
		#plt.setp( axes[0].get_yticklabels(), visible=False)

		axes[1].plot(star_track[:,0], res_list_d)
		axes[1].set_ylabel('Response')
		axes[1].set_xlabel('time/s')

		plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
		            wspace=0, hspace=0)
		plt.savefig('../plots/%s/star/res/2p/res%d-%d_track.png'%(scan_name,star_num,order),dpi=190)  
		plt.clf()


def makeGaussian(size, fwhm = 1, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

def split_into_chunks(asp, tranges):
	asp_t = (asp[:,0]*1000).astype(int)


def downsample(ar, fact):
	sx, sy = ar.shape
	X, Y = np.ogrid[0:sx, 0:sy]
	regions = sy/fact * (X/fact) + Y/fact
	res = ndimage.sum(ar, labels=regions, index=np.arange(regions.max() + 1))
	res.shape = (sx/fact, sy/fact)
	return res

def run_model_chunk(pid, scan_name, chunk, cata_target, cata_lists, size_list, cata_target_len, dead_tmp, asp_uni, return_dict):
	count_list=[]
	exp_list=[]
	chunk_num = 0
	for trange in chunk:
		output = '../fits/%s/extra/new_black/exposure_chunk%d.fits'%(scan_name, chunk_num)
		if  os.path.exists(output):
			print 'exists %d'%(chunk_num)
			chunk_num+=1
			continue
		t_mask = np.zeros(asp_uni.shape[0])
		begin = np.where((asp_uni[:,0]*1000).astype(int) - trange[0] == 0)[0] 
		end = np.where((asp_uni[:,0]*1000).astype(int) - trange[1] == 0)[0]
		t_mask[begin:end+1] = 1

		t_asp = asp_uni[t_mask>0]
		t_dead = dead_tmp[t_mask>0]

		t_length = end+1-begin
		print('time length:%d, %f, %f'%(t_length, trange[0], trange[1]))
		count = np.zeros([t_length,800,800])
		exp_count = np.zeros([t_length,800,800])
		for star_num in range(cata_target_len):
			#print star_num
			star_co = cata_target[star_num, 0:2]
			star_flux = cata_target[star_num, 3]
			#print star_flux

			ra = np.ones(t_asp.shape[0])*star_co[0]
			dec = np.ones(t_asp.shape[0])*star_co[1]

			xi, eta = gn.gnomfwd_simple(ra, dec, t_asp[:,1], t_asp[:,2], -t_asp[:,3], 1/36000., 0.0)

			coo = np.array([xi, eta]).T
			coo = ((coo/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800

			for i in [-1,0,1]:
				for j in [-1,0,1]:
					coo_tmp = np.zeros([t_length, 3])
					coo_tmp[:,0] = np.arange(t_length)
					coo_tmp[:,1] = coo[:,0]+i
					coo_tmp[:,2] = coo[:,1]+j
					
					mask_0 = np.logical_and(coo_tmp[:,1]>=0, coo_tmp[:,1]<800)
					mask_1 = np.logical_and(coo_tmp[:,2]>=0, coo_tmp[:,2]<800)
					mask = np.logical_and(mask_0, mask_1)
					
					total = 1+4*0.0625+4*0.00390625
					coo_tmp = coo_tmp[mask]
					t_dead_tmp = t_dead[mask]
					coo_tmp = np.array(coo_tmp, dtype=int)
					factor = 1.
					if (i,j) in set([(-1,0),(1,0),(0,-1),(0,1)]):
						factor = 0.0625
					elif (i,j) in set([(-1,-1),(1,-1),(1,1),(-1,1)]):
						factor = 0.00390625
					else:
						factor = 1.
					#count[coo_tmp[:,0], coo_tmp[:,1]] += 0.005*factor*star_flux/total
					count[coo_tmp[:,0], coo_tmp[:,1], coo_tmp[:,2]] += factor*star_flux/total

		for mask_num in range(3):
			mask_cata = cata_lists[mask_num]
			size = size_list[mask_num]
			radius = size[-1]+0.5
			for star_num in range(mask_cata.shape[0]):
				#print star_num
				star_co = mask_cata[star_num, 0:2]
				#print star_flux

				ra = np.ones(t_asp.shape[0])*star_co[0]
				dec = np.ones(t_asp.shape[0])*star_co[1]

				xi, eta = gn.gnomfwd_simple(ra, dec, t_asp[:,1], t_asp[:,2], -t_asp[:,3], 1/36000., 0.0)

				coo = np.array([xi, eta]).T
				coo = ((coo/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800

				for i in size:
					for j in size:
						if i**2+j**2>radius**2:
							continue
						coo_tmp = np.zeros([t_length, 3])
						coo_tmp[:,0] = np.arange(t_length)
						coo_tmp[:,1] = coo[:,0]+i
						coo_tmp[:,2] = coo[:,1]+j
						
						mask_0 = np.logical_and(coo_tmp[:,1]>=0, coo_tmp[:,1]<800)
						mask_1 = np.logical_and(coo_tmp[:,2]>=0, coo_tmp[:,2]<800)
						mask = np.logical_and(mask_0, mask_1)
						coo_tmp = coo_tmp[mask]

						coo_tmp = np.array(coo_tmp, dtype=int)

						exp_count[coo_tmp[:,0], coo_tmp[:,1], coo_tmp[:,2]] -= 1
		for t_step in range(t_length): 
			exp_count[t_step][exp_count[t_step]>=0] = 0.005*(1-t_dead[t_step])
		exp_count[exp_count<0] = 0.
		count *= exp_count

		exp_count_sum = np.sum(exp_count, axis=0)
		count_sum = np.sum(count, axis=0)
		exp_list.append(exp_count_sum)
		count_list.append(count_sum)
		exp_count=count=None
		gc.collect()
		chunk_num += 1
	return_dict[pid] = (count_list, exp_list)



if __name__ == '__main__':
	if False:	
		scan_name = 'AIS_GAL_SCAN_00014_0001'
		manager = Manager()
		return_dict = manager.dict()
		num_list = [0, 1, 2, 3, 4, 5, 6, 9, 13, 39, 44, 49, 56, 57, 60, 75, 85, 88, 100, 118, 126]
		for num in num_list:
			print num
			res_star_new(0, scan_name, num, num, return_dict)

	if False:
		scan_name = sys.argv[1]
		#scan_name = 'AIS_GAL_SCAN_00014_0001'

		output = '../fits/%s/extra/flat.fits'%(scan_name)
		print output
		path = os.path.dirname(output)
		if not os.path.exists(path):
			os.makedirs(path)
		else:
			print 'exists'

		cata = spi.load_obj('../data/%s_starset_extra_all_star'%scan_name)
		cata_a = np.array(list(cata))
		cata_a = cata_a[np.logical_and(cata_a[:,4]>15.5,cata_a[:,4]<17.5)]
		cata_len = len(cata_a)
		cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)

		print cata_len

		asp = np.load('../data/photon_list/%s_asp_cal.npy'%scan_name)
		asp_uni, uni_index=np.unique(asp[:,0], True)
		asp_uni = asp[uni_index]
		ix_tmp = (np.round(asp_uni[:,0]-asp_uni[0,0])+1).astype(int)
		dead_tmp = np.zeros(ix_tmp.shape[0])
		scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%scan_name)[1].data
		limit = scst_tmp['t_dead_nuv'].shape[0]-1
		ix_tmp[ix_tmp>limit] = limit
		dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]

		count = np.zeros([800,800])
		#for star_num in xrange(1,10):
		for star_num in range(cata_len):
			print star_num
			star_co = cata_a[star_num, 0:2]
			star_flux = cata_a[star_num, 3]
			print star_flux

			ra = np.ones(asp_uni.shape[0])*star_co[0]
			dec = np.ones(asp_uni.shape[0])*star_co[1]

			xi, eta = gn.gnomfwd_simple(ra, dec, asp_uni[:,1], asp_uni[:,2], -asp_uni[:,3], 1/36000., 0.0)

			coo = np.array([xi, eta]).T
			#print coo.shape
			#print coo
			coo = ((coo/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800

			'''
			mask_0 = np.logical_and(coo[:,0]>=0, coo[:,0]<800)
			mask_1 = np.logical_and(coo[:,1]>=0, coo[:,1]<800)

			mask = np.logical_and(mask_0, mask_1)

			coo = coo[mask]
			dead_t = dead_tmp[mask]

			flat[coo[:,0], coo[:,1]] += 0.005*(1-dead_t)
			'''

			for i in [-1,0,1]:
				for j in [-1,0,1]:
					coo_tmp = np.zeros(coo.shape)
					coo_tmp[:,0] = coo[:,0]+i
					coo_tmp[:,1] = coo[:,1]+j
					mask = (coo_tmp[:,0]-400)**2+(coo_tmp[:,1]-400)**2<148996
					'''
					mask_0 = np.logical_and(coo_tmp[:,0]>=0, coo_tmp[:,0]<800)
					mask_1 = np.logical_and(coo_tmp[:,1]>=0, coo_tmp[:,1]<800)
					mask = np.logical_and(mask_0, mask_1)
					'''
					coo_tmp = coo_tmp[mask]
					dead_t = dead_tmp[mask]
					coo_tmp = np.array(coo_tmp, dtype=int)
					if i==0 and j==0:
						count[coo_tmp[:,0], coo_tmp[:,1]] += 0.005*(1-dead_t)*star_flux
					else:
						count[coo_tmp[:,0], coo_tmp[:,1]] += 0.005*(1-dead_t)/2.*star_flux

		flat = pyfits.open('../data/cal/NUV_flat.fits')
		hdu = pyfits.PrimaryHDU(count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/%s/extra/star.fits'%(scan_name), clobber=False)

#flat gaussian
	if False:
		scan_name = sys.argv[1]
		#scan_name = 'AIS_GAL_SCAN_00014_0001'

		output = '../fits/%s/extra/flat.fits'%(scan_name)
		print output
		path = os.path.dirname(output)
		if not os.path.exists(path):
			os.makedirs(path)
		else:
			print 'exists'

		cata = spi.load_obj('../data/%s_starset_extra_all_star'%scan_name)
		cata_a = np.array(list(cata))
		cata_a = cata_a[np.logical_and(cata_a[:,4]>15.5,cata_a[:,4]<17.5)]
		cata_len = len(cata_a)
		cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)

		print cata_len

		asp = np.load('../data/photon_list/%s_asp_cal.npy'%scan_name)
		asp_uni, uni_index=np.unique(asp[:,0], True)
		asp_uni = asp[uni_index]
		ix_tmp = (np.round(asp_uni[:,0]-asp_uni[0,0])+1).astype(int)
		dead_tmp = np.zeros(ix_tmp.shape[0])
		scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%scan_name)[1].data
		limit = scst_tmp['t_dead_nuv'].shape[0]-1
		ix_tmp[ix_tmp>limit] = limit
		dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]

		count = np.zeros([800,800])
		#for star_num in xrange(1,10):
		for star_num in range(cata_len):
			print star_num
			star_co = cata_a[star_num, 0:2]
			star_flux = cata_a[star_num, 3]
			print star_flux

			ra = np.ones(asp_uni.shape[0])*star_co[0]
			dec = np.ones(asp_uni.shape[0])*star_co[1]

			xi, eta = gn.gnomfwd_simple(ra, dec, asp_uni[:,1], asp_uni[:,2], -asp_uni[:,3], 1/36000., 0.0)

			coo = np.array([xi, eta]).T
			#print coo.shape
			#print coo
			coo = ((coo/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800

			'''
			mask_0 = np.logical_and(coo[:,0]>=0, coo[:,0]<800)
			mask_1 = np.logical_and(coo[:,1]>=0, coo[:,1]<800)

			mask = np.logical_and(mask_0, mask_1)

			coo = coo[mask]
			dead_t = dead_tmp[mask]

			flat[coo[:,0], coo[:,1]] += 0.005*(1-dead_t)
			'''

			for i in [-2,-1,0,1,2]:
				for j in [-2,-1,0,1,2]:
					coo_tmp = np.zeros(coo.shape)
					coo_tmp[:,0] = coo[:,0]+i
					coo_tmp[:,1] = coo[:,1]+j
					mask = (coo_tmp[:,0]-400)**2+(coo_tmp[:,1]-400)**2<148996
					'''
					mask_0 = np.logical_and(coo_tmp[:,0]>=0, coo_tmp[:,0]<800)
					mask_1 = np.logical_and(coo_tmp[:,1]>=0, coo_tmp[:,1]<800)
					mask = np.logical_and(mask_0, mask_1)
					'''
					coo_tmp = coo_tmp[mask]
					dead_t = dead_tmp[mask]
					coo_tmp = np.array(coo_tmp, dtype=int)
					factor = 1.
					if (i,j) in set([(-1,0),(1,0),(0,-1),(0,1)]):
						factor = 0.73486724613779941
					elif (i,j) in set([(-1,-1),(1,-1),(1,1),(-1,1)]):
						factor = 0.54002986944615305
					elif (i,j) in set([(-2,-2),(2,-2),(2,2),(-2,2)]):
						factor = 0.085049375010898556
					elif (i,j) in set([(-2,-1),(-2,1),(2,-1),(2,1),(-1,-2),(1,-2),(-1,2),(1,2)]):
						factor = 0.21431099571326823
					elif (i,j) in set([(-2,0),(2,0),(0,-2),(0,2)]):
						factor = 0.29163225989402913
					else:
						factor = 1.
					count[coo_tmp[:,0], coo_tmp[:,1]] += 0.005*(1-dead_t)*star_flux*factor

		flat = pyfits.open('../data/cal/NUV_flat.fits')
		hdu = pyfits.PrimaryHDU(count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/%s/extra/star_gaussian.fits'%(scan_name), clobber=False)

	if False:
		scan_name = 'AIS_GAL_SCAN_00014_0001'
		flat = pyfits.open('../fits/%s/extra/flat.fits'%(scan_name))
		star = pyfits.open('../fits/%s/extra/star_gaussian.fits'%(scan_name))
		star_data = star[0].data
		star_data[np.where(star_data<=0.00001)] = 1.
		star_data = np.swapaxes(star_data,0,1)
		count = flat[0].data/star_data

		hdu = pyfits.PrimaryHDU(count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/%s/extra/NUV_flat_gaussian.fits'%(scan_name), clobber=False)

	if False:
		name_list = ['AIS_GAL_SCAN_00023_0001', 'AIS_GAL_SCAN_00032_0001', 'AIS_GAL_SCAN_00041_0001','AIS_GAL_SCAN_00050_0001']
		scan_name = 'AIS_GAL_SCAN_00050_0001'
		count = np.zeros([800,800])
		star_count = np.zeros([800,800])
		for scan_name in name_list:
			flat = pyfits.open('../fits/%s/extra/flat.fits'%(scan_name))
			star = pyfits.open('../fits/%s/extra/star.fits'%(scan_name))
			star_data = star[0].data
			count += flat[0].data
			star_count += star_data
		star_count[np.where(star_count<=0.00001)] = 1.
		'''
		star_count = np.swapaxes(star_count,0,1)
		count = count/star_count
		hdu = pyfits.PrimaryHDU(count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/flat/NUV_flat_no_14.fits', clobber=False)
		'''
		hdu = pyfits.PrimaryHDU(count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/flat/flat_no_14.fits', clobber=False)

		hdu = pyfits.PrimaryHDU(star_count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/flat/star_no_14.fits', clobber=False)

#flat gaussian chunks
	if False:
		scan_name = sys.argv[1]
		#scan_name = 'AIS_GAL_SCAN_00014_0001'

		output = '../fits/%s/extra/new_faint/flat.fits'%(scan_name)
		print output
		path = os.path.dirname(output)
		if not os.path.exists(path):
			os.makedirs(path)
		else:
			print 'exists'

		cata = spi.load_obj('../data/%s_starset_extra_all_star'%scan_name)
		cata_a = np.array(list(cata))
		cata_a = cata_a[np.logical_and(cata_a[:,4]>16.0,cata_a[:,4]<18.0)]
		cata_len = len(cata_a)
		cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)

		tranges = np.load('../fits/%s/extra/new_faint/tranges.npy'%scan_name)

		print cata_len

		asp = np.load('../data/photon_list/%s_asp_cal.npy'%scan_name)
		asp_uni, uni_index=np.unique(asp[:,0], True)
		asp_uni = asp[uni_index]
		ix_tmp = (np.round(asp_uni[:,0]-asp_uni[0,0])+1).astype(int)
		dead_tmp = np.zeros(ix_tmp.shape[0])
		scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%scan_name)[1].data
		limit = scst_tmp['t_dead_nuv'].shape[0]-1
		ix_tmp[ix_tmp>limit] = limit
		dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]

		flat = pyfits.open('../data/cal/NUV_flat.fits')
		chunk_num = 0
		for trange in tranges:
			output = '../fits/%s/extra/new_faint/exposure_chunk%d.fits'%(scan_name, chunk_num)
			if  os.path.exists(output):
				print 'exists %d'%(chunk_num)
				chunk_num+=1
				continue
			t_mask = np.zeros(asp_uni.shape[0])
			begin = np.where((asp_uni[:,0]*1000).astype(int) - trange[0] == 0)[0] 
			end = np.where((asp_uni[:,0]*1000).astype(int) - trange[1] == 0)[0]
			t_mask[begin:end+1] = 1

			t_asp = asp_uni[t_mask>0]
			t_dead = dead_tmp[t_mask>0]

			count = np.zeros([800,800])
			exp_count = np.zeros([800,800])
			#for star_num in xrange(1,10):
			for star_num in range(cata_len):
				#print star_num
				star_co = cata_a[star_num, 0:2]
				star_flux = cata_a[star_num, 3]
				#print star_flux

				ra = np.ones(t_asp.shape[0])*star_co[0]
				dec = np.ones(t_asp.shape[0])*star_co[1]

				xi, eta = gn.gnomfwd_simple(ra, dec, t_asp[:,1], t_asp[:,2], -t_asp[:,3], 1/36000., 0.0)

				coo = np.array([xi, eta]).T
				#print coo.shape
				#print coo
				coo = ((coo/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800

				for i in [-1,0,1]:
					for j in [-1,0,1]:
						coo_tmp = np.zeros(coo.shape)
						coo_tmp[:,0] = coo[:,0]+i
						coo_tmp[:,1] = coo[:,1]+j
						
						mask_0 = np.logical_and(coo_tmp[:,0]>=0, coo_tmp[:,0]<800)
						mask_1 = np.logical_and(coo_tmp[:,1]>=0, coo_tmp[:,1]<800)
						mask = np.logical_and(mask_0, mask_1)
						
						total = 1+4*0.8133*0.017869487123182051+4*0.3322*0.0003193185700455692
						coo_tmp = coo_tmp[mask]
						t_dead_tmp = t_dead[mask]
						coo_tmp = np.array(coo_tmp, dtype=int)
						factor = 1.
						if (i,j) in set([(-1,0),(1,0),(0,-1),(0,1)]):
							factor = 0.8133*0.017869487123182051
						elif (i,j) in set([(-1,-1),(1,-1),(1,1),(-1,1)]):
							factor = 0.3322*0.0003193185700455692
						else:
							factor = 1.
						#count[coo_tmp[:,0], coo_tmp[:,1]] += 0.005*factor*star_flux/total
						count[coo_tmp[:,0], coo_tmp[:,1]] += 0.005*(1-t_dead_tmp)*factor*star_flux/total
						exp_count[coo_tmp[:,0], coo_tmp[:,1]] += 0.005*(1-t_dead_tmp)

			hdu = pyfits.PrimaryHDU(count)
			hdu.header = flat[0].header
			hdulist = pyfits.HDUList([hdu])
			hdulist.writeto('../fits/%s/extra/new_faint/star_chunk%d.fits'%(scan_name, chunk_num), clobber=False)

			hdu = pyfits.PrimaryHDU(exp_count)
			hdu.header = flat[0].header
			hdulist = pyfits.HDUList([hdu])
			hdulist.writeto('../fits/%s/extra/new_faint/exposure_chunk%d.fits'%(scan_name, chunk_num), clobber=False)
			chunk_num += 1


#flat gaussian chunks, black
	if True:
		scan_name = sys.argv[1]
		#scan_name = 'AIS_GAL_SCAN_00014_0001'

		output = '../fits/%s/extra/new_black/flat.fits'%(scan_name)
		print output
		path = os.path.dirname(output)
		if not os.path.exists(path):
			os.makedirs(path)
		else:
			print 'exists'

		cata = spi.load_obj('../data/%s_starset_extra_all_star'%scan_name)
		cata_a = np.array(list(cata))
		mask_faint = np.logical_and(cata_a[:,3]>6, cata_a[:,3]<40)
		mask_mid = np.logical_and(cata_a[:,3]>=40, cata_a[:,3]<120)
		mask_bright = np.logical_and(cata_a[:,3]>=120, cata_a[:,3]<250)
		mask_ultra = cata_a[:,3]>=250
		mask_list = [mask_mid, mask_bright, mask_ultra]
		size_list = [[-3,-2,-1,0,1,2,3],[-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9], 
		            [-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]]
		cata_len = len(cata_a)
		cata_target = cata_a[mask_faint]
		cata_lists = []
		for mask in mask_list:
			cata_lists.append(cata_a[mask])

		tranges = np.load('../fits/%s/extra/new_black/tranges.npy'%scan_name)

		cata_target_len = cata_target.shape[0]
		print cata_target_len

		asp = np.load('../data/photon_list/%s_asp_cal_inter_half_sec.npy'%scan_name)
		asp_uni, uni_index=np.unique(asp[:,0], True)
		asp_uni = asp[uni_index]
		ix_tmp = (np.round(asp_uni[:,0]-asp_uni[0,0])+1).astype(int)
		dead_tmp = np.zeros(ix_tmp.shape[0])
		scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%scan_name)[1].data
		limit = scst_tmp['t_dead_nuv'].shape[0]-1
		ix_tmp[ix_tmp>limit] = limit
		dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]

		flat = pyfits.open('../data/cal/NUV_flat.fits')
		count_list=[]
		exp_list=[]

		#total_p_num = 10
		total_p_num = int(sys.argv[2])

		manager = Manager()
		return_dict = manager.dict()
		p_list = []

		pid = 0
		chunks = split_seq(tranges, total_p_num)
		print len(chunks)
		for chunk in chunks:
			print chunk.shape
			p = Process(target=run_model_chunk, args=(pid, scan_name, chunk, cata_target, cata_lists, size_list, cata_target_len, dead_tmp, asp_uni, return_dict))
			p.start()
			p_list.append(p)
			pid += 1
		print pid

		for p in p_list:
			p.join()
		print 'all done'
		for i in range(0, pid):
			counts, exps = return_dict[i]
			count_list += counts
			exp_list += exps
		np.save('../fits/%s/extra/new_black/exposure.fits'%(scan_name), exp_list)
		np.save('../fits/%s/extra/new_black/star.fits'%(scan_name), count_list)


#flat exposure chunks
	if False:
		scan_name = sys.argv[1]
		#scan_name = 'AIS_GAL_SCAN_00014_0001'

		output = '../fits/%s/extra/new_1280/flat.fits'%(scan_name)
		print output
		path = os.path.dirname(output)
		if not os.path.exists(path):
			os.makedirs(path)
		else:
			print 'exists'

		cata = spi.load_obj('../data/%s_starset_extra_all_star'%scan_name)
		cata_a = np.array(list(cata))
		cata_a = cata_a[np.logical_and(cata_a[:,4]>16.0,cata_a[:,4]<18.0)]
		cata_len = len(cata_a)
		cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)

		tranges = np.load('../fits/%s/extra/new_1280/tranges.npy'%scan_name)

		print cata_len

		asp = np.load('../data/photon_list/%s_asp_cal.npy'%scan_name)
		asp_uni, uni_index=np.unique(asp[:,0], True)
		asp_uni = asp[uni_index]
		ix_tmp = (np.round(asp_uni[:,0]-asp_uni[0,0])+1).astype(int)
		dead_tmp = np.zeros(ix_tmp.shape[0])
		scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%scan_name)[1].data
		limit = scst_tmp['t_dead_nuv'].shape[0]-1
		ix_tmp[ix_tmp>limit] = limit
		dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]

		chunk_num = 0
		for trange in tranges:
			t_mask = np.zeros(asp_uni.shape[0])
			begin = np.where((asp_uni[:,0]*1000).astype(int) - trange[0] == 0)[0] 
			end = np.where((asp_uni[:,0]*1000).astype(int) - trange[1] == 0)[0]
			t_mask[begin:end+1] = 1

			t_asp = asp_uni[t_mask>0]
			t_dead = dead_tmp[t_mask>0]

			count = np.zeros([800,800])
			#for star_num in xrange(1,10):
			for star_num in range(cata_len):
				print star_num
				star_co = cata_a[star_num, 0:2]
				star_flux = cata_a[star_num, 3]
				print star_flux

				ra = np.ones(t_asp.shape[0])*star_co[0]
				dec = np.ones(t_asp.shape[0])*star_co[1]

				xi, eta = gn.gnomfwd_simple(ra, dec, t_asp[:,1], t_asp[:,2], -t_asp[:,3], 1/36000., 0.0)

				coo = np.array([xi, eta]).T
				#print coo.shape
				#print coo
				coo = ((coo/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*800

				for i in [-1,0,1]:
					for j in [-1,0,1]:
						coo_tmp = np.zeros(coo.shape)
						coo_tmp[:,0] = coo[:,0]+i
						coo_tmp[:,1] = coo[:,1]+j
						
						mask_0 = np.logical_and(coo_tmp[:,0]>=0, coo_tmp[:,0]<800)
						mask_1 = np.logical_and(coo_tmp[:,1]>=0, coo_tmp[:,1]<800)
						mask = np.logical_and(mask_0, mask_1)
						
						coo_tmp = coo_tmp[mask]
						t_dead_tmp = t_dead[mask]
						coo_tmp = np.array(coo_tmp, dtype=int)
						factor = 1.
						count[coo_tmp[:,0], coo_tmp[:,1]] += 0.005*(1-t_dead_tmp)*factor

			flat = pyfits.open('../data/cal/NUV_flat.fits')
			hdu = pyfits.PrimaryHDU(count)
			hdu.header = flat[0].header
			hdulist = pyfits.HDUList([hdu])
			hdulist.writeto('../fits/%s/extra/new_1280/exposure_chunk%d.fits'%(scan_name, chunk_num), clobber=False)
			chunk_num += 1

