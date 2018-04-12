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

if __name__ == '__main__':

#Fit flat downsample, black
	if False:
		scan_list = ['AIS_GAL_SCAN_02057_0002', 'AIS_GAL_SCAN_02057_0003', 'AIS_GAL_SCAN_02066_0003', 'AIS_GAL_SCAN_02075_0002']
		#scan_list = ['AIS_GAL_SCAN_00176_0001', 'AIS_GAL_SCAN_00185_0001', 'AIS_GAL_SCAN_00194_0001', 'AIS_GAL_SCAN_00203_0001']
		
		list_index = int(sys.argv[1])
		size = int(sys.argv[2])
		print list_index
		with open('../fits/flat/chunk_list/%d'%list_index) as f:
			scan_list = f.read().splitlines()
		print scan_list
		
		dirname = 'black_f8_%d'%list_index
		print dirname
		flat_o = pyfits.open('../data/cal/NUV_flat.fits')
		flat_list = []
		star_list = []
		exposure_list = []
		scan_num = 0
		#size = 100
		chunk_size = 2560*100*100/size/size#2560#5120*25*25/size/size
		part_step = 5120/chunk_size


		output = '../fits/flat/chunk/%s/downsample_%d/NUV_flat_full_%d_new.npy'%(dirname, size, size)
		print output
		path = os.path.dirname(output)
		if not os.path.exists(path):
			os.makedirs(path)
		else:
			print 'exists'

		for scan_name in scan_list:
			print scan_name
			data = np.load('../fits/%s/extra/new_black/count_list.npy'%(scan_name))
			star = np.load('../fits/%s/extra/new_black/star.fits.npy'%(scan_name))
			exp = np.load('../fits/%s/extra/new_black/exposure.fits.npy'%(scan_name))

			for part in range(chunk_size/4):
				flat_tmp = np.zeros([800,800])
				star_tmp = np.zeros([800,800])
				exposure_tmp = np.zeros([800,800])
				for chunk_num in range(part*part_step, part*part_step+part_step):
					print chunk_num
					flat_tmp += np.swapaxes(data[chunk_num],0,1)

					star_tmp += star[chunk_num]

					exposure_tmp +=exp[chunk_num]

				down_fact = 800/size
				flat_tmp = downsample(flat_tmp, down_fact)
				star_tmp = downsample(star_tmp, down_fact)
				exposure_tmp = downsample(exposure_tmp, down_fact)

				flat_list.append(flat_tmp.flatten())
				star_list.append(star_tmp.flatten())
				exposure_list.append(exposure_tmp.flatten())

			scan_num += 1
			flat_tmp=star_tmp=exposure_tmp=data=star=exp=None
			gc.collect()

		np.save('../fits/flat/chunk/%s/downsample_%d/data_chunk_%d.npy'%(dirname, size, size), flat_list)
		np.save('../fits/flat/chunk/%s/downsample_%d/star_chunk_%d.npy'%(dirname, size, size), star_list)
		np.save('../fits/flat/chunk/%s/downsample_%d/exp_chunk_%d.npy'%(dirname, size, size), exposure_list)

		print('data chunked')


		field = np.zeros([size,size])
		field = field.flatten()

		flat = np.array(flat_list)
		star = np.array(star_list)
		bkg = np.array(exposure_list)
		flat_list=star_list=exposure_list=None
		gc.collect()

		coefficients = []
		print flat.shape
		print star.shape
		print bkg.shape
		count = np.zeros([size,size])
		#star_count = np.zeros([size,size])
		#bkg_count = np.zeros([size,size])
		model = np.zeros([chunk_size, size, size])
		bkg_flat = np.zeros([size,size])
		for i in range(0, size):
			for j in range(0, size):
				if (i-size/2.)**2+(j-size/2.)**2>148996*(size/800.)**2:
					continue
				field[i*size+j] = 1
				y = flat[:,field>0]
				s = star[:,field>0]
				b = bkg[:,field>0]
				x = np.concatenate([s, b], axis=1)

				result = lss.linear_least_squares(x, y)
				count[i,j] = result[0]
				bkg_flat[i,j] = result[1]
				coefficients.append(result)

				model[:,i,j] = (s*result[0])[:,0]+(b*result[1])[:,0]
				#star_count[i,j] = np.sum((s*result[0])[:,0], axis=0)
				#bkg_count[i,j] = np.sum((b*result[1])[:,0], axis=0)

				field[i*size+j] = 0

		coefficients = np.array(coefficients)
		np.save('../fits/flat/chunk/%s/downsample_%d/NUV_flat_full_%d_new.npy'%(dirname, size, size), coefficients)

		flat = pyfits.open('../data/cal/NUV_flat.fits')
		
		hdu = pyfits.PrimaryHDU(model)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/flat/chunk/%s/downsample_%d/model_count.fits'%(dirname, size), clobber=False)


		hdu = pyfits.PrimaryHDU(count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/flat/chunk/%s/downsample_%d/NUV_flat_try_full_%d_nocon_new.fits'%(dirname, size,size), clobber=False)

		'''
		hdu = pyfits.PrimaryHDU(star_count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/flat/chunk/%s/downsample_%d/star_count.fits'%(dirname, size), clobber=False)

		hdu = pyfits.PrimaryHDU(bkg_count)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/flat/chunk/%s/downsample_%d/bkg_count.fits'%(dirname, size), clobber=False)
		'''

		hdu = pyfits.PrimaryHDU(bkg_flat)
		hdu.header = flat[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/flat/chunk/%s/downsample_%d/bkg_flat.fits'%(dirname, size), clobber=False)


		'''
		for i in range(chunk_size):
			hdu = pyfits.PrimaryHDU(star_count[i])
			hdu.header = flat[0].header
			hdulist = pyfits.HDUList([hdu])
			hdulist.writeto('../fits/flat/chunk/downsample_%d/model_count%d.fits'%(size, i), clobber=False)
		'''


#Fit flat downsample, black, plot
	if False:
		list_index = int(sys.argv[1])
		dirname = 'black_f8_2'
		dirname = 'black_f8_%d'%list_index
		print dirname
		
		for size in [ 100]:
			#star = pyfits.open('../fits/flat/chunk/%s/downsample_%d/star_count.fits'%(dirname, size))[0].data
			#bkg = pyfits.open('../fits/flat/chunk/%s/downsample_%d/bkg_count.fits'%(dirname, size))[0].data
			#total = star+bkg
			bkg_flat = pyfits.open('../fits/flat/chunk/%s/downsample_%d/bkg_flat.fits'%(dirname, size))[0].data
			star_flat = pyfits.open('../fits/flat/chunk/%s/downsample_%d/NUV_flat_try_full_%d_nocon_new.fits'%(dirname, size, size))[0].data

			'''
			imgplot=plt.imshow(star, origin='lower', interpolation='none')
			imgplot.set_cmap('gray')
			plt.title('star_model')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/star_model_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

			imgplot=plt.imshow(bkg, origin='lower', interpolation='none')
			imgplot.set_cmap('gray')
			plt.title('bg model')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/bg_model_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

			imgplot=plt.imshow(total, origin='lower', interpolation='none')
			imgplot.set_cmap('gray')
			plt.title('model')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/model_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()
			'''

			imgplot=plt.imshow(bkg_flat, origin='lower', interpolation='none')
			imgplot.set_cmap('gray')
			plt.title('bg_flat X bg_flux')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/bg_flat_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

			imgplot=plt.imshow(star_flat, origin='lower', interpolation='none')
			imgplot.set_cmap('gray')
			plt.title('star_flat')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/star_flat_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()


#Fit flat downsample, black, chunk plot
	if False:
		#dirname = 'black_f8_2'
		list_index = int(sys.argv[1])
		dirname = 'black_f8_%d'%list_index
		print dirname
		
		for size in [ 100]:
			chunk_num = 160*(400/size)**2
			total = pyfits.open('../fits/flat/chunk/%s/downsample_%d/model_count.fits'%(dirname, size))[0].data
			data = np.load('../fits/flat/chunk/%s/downsample_%d/data_chunk_%d.npy'%(dirname, size, size))
			star = np.load('../fits/flat/chunk/%s/downsample_%d/star_chunk_%d.npy'%(dirname, size, size))
			exp = np.load('../fits/flat/chunk/%s/downsample_%d/exp_chunk_%d.npy'%(dirname, size, size))

			print data.shape
			print total.shape
			imgplot=plt.imshow(data[chunk_num/1.5].reshape([size, size]), origin='lower', interpolation='none', vmin=0, vmax=40)
			imgplot.set_cmap('gray')
			plt.title('Data')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/_data_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

			imgplot=plt.imshow(star[chunk_num/1.5].reshape([size, size]), origin='lower', interpolation='none')
			imgplot.set_cmap('gray')
			plt.title('star_flux X exposure_time')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/_star_flux_exp_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

			imgplot=plt.imshow(exp[chunk_num/1.5].reshape([size, size]), origin='lower', interpolation='none')
			imgplot.set_cmap('gray')
			plt.title('exposure_time')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/_exposure_time_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

			imgplot=plt.imshow(total[chunk_num/1.5], origin='lower', interpolation='none', vmin=0, vmax=40)
			imgplot.set_cmap('gray')
			plt.title('model')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/_model_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

#compare black, plot
	if False:
		dirnames = ['black_f8', 'black_f8_1']
		
		star_flats=[]
		bkg_flats=[]
		for size in [400, 200, 100]:
			star_flats=[]
			bkg_flats=[]
			for dirname in dirnames:
				bkg_flats.append(pyfits.open('../fits/flat/chunk/%s/downsample_%d/bkg_flat.fits'%(dirname, size))[0].data)
				star_flats.append(pyfits.open('../fits/flat/chunk/%s/downsample_%d/NUV_flat_try_full_%d_nocon_new.fits'%(dirname, size, size))[0].data)

			star_ratio = star_flats[0]/star_flats[1]
			bkg_ratio = bkg_flats[0]/bkg_flats[1]

			imgplot=plt.imshow(bkg_ratio, origin='lower', interpolation='none', vmin=1.5, vmax=2.2)
			imgplot.set_cmap('gray')
			plt.title('bg ratio')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/bg_ratio_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

			imgplot=plt.imshow(star_ratio, origin='lower', interpolation='none', vmin=0.9, vmax=1.1)
			imgplot.set_cmap('gray')
			plt.title('star ratio')
			plt.colorbar()
			plt.tight_layout()
			plt.savefig('../fits/flat/chunk/%s/downsample_%d/star_ratio_%d.png'%(dirname, size, size), dpi=190)
			plt.clf()

	if False:
		star_flats = []
		bkg_flats = [] 
		size=100
		for i in range(9):
			dirname = 'black_f8_%d'%i
			star_flats.append(pyfits.open('../fits/flat/chunk/%s/downsample_%d/NUV_flat_try_full_%d_nocon_new.fits'%(dirname, size, size))[0].data)
			bkg_flats.append(pyfits.open('../fits/flat/chunk/%s/downsample_%d/bkg_flat.fits'%(dirname, size))[0].data)

		m_flat = np.median(star_flats, axis=0)
		print m_flat.shape
		dev = np.square(star_flats-m_flat)
		m_dev = np.mean(dev, axis=0)


		m_flat_bkg = np.median(bkg_flats, axis=0)
		dev_bkg = np.square(bkg_flats-m_flat_bkg)
		m_dev_bkg = np.mean(dev_bkg, axis=0)

		plt.imshow(m_dev, origin='lower', interpolation='none')
		plt.set_cmap('gray')
		plt.colorbar()
		plt.savefig('../plots/flat/var_map100.png', dpi=190)
		plt.show()

		plt.imshow(m_flat, origin='lower', interpolation='none')
		plt.set_cmap('gray')
		plt.colorbar()
		plt.savefig('../plots/flat/median_flat100.png', dpi=190)
		plt.show()

		plt.imshow(m_dev_bkg, origin='lower', interpolation='none')
		plt.set_cmap('gray')
		plt.colorbar()
		plt.savefig('../plots/flat/var_map_bkg100.png', dpi=190)
		plt.show()

		plt.imshow(m_flat_bkg, origin='lower', interpolation='none')
		plt.set_cmap('gray')
		plt.colorbar()
		plt.savefig('../plots/flat/median_bkg_flat100.png', dpi=190)
		plt.show()


	if True:
		star_flats = []
		bkg_flats = [] 
		star_medians = []
		bkg_medians = []
		size=100
		margin = 25
		for i in range(9):
			dirname = 'black_f8_%d'%i
			star_flats.append(pyfits.open('../fits/flat/chunk/%s/downsample_%d/NUV_flat_try_full_%d_nocon_new.fits'%(dirname, size, size))[0].data)
			bkg_flats.append(pyfits.open('../fits/flat/chunk/%s/downsample_%d/bkg_flat.fits'%(dirname, size))[0].data)
			star_medians.append(np.median(star_flats[i][margin:size-margin, margin:size-margin]))
			bkg_medians.append(np.median(bkg_flats[i][margin:size-margin, margin:size-margin]))

		star_median = np.median(star_medians)
		bkg_median = 1.#np.median(bkg_medians)
		print star_median
		print bkg_median
		for i in range(9):
			star_flats[i] = star_flats[i]*star_median/star_medians[i]
			bkg_flats[i] = bkg_flats[i]*bkg_median/bkg_medians[i]
		
		flat_max = np.max(star_flats)

		'''
		flat_max = np.max(star_flats)
		f, ax = plt.subplots(9,9)
		f.set_size_inches(12.9,12.9)
		for i in range(9):
			for j in range(i+1):
				ax[i,j].plot(star_flats[i].flatten(), star_flats[j].flatten(), '.k')
				ax[i,j].set_xlim(-0.1, flat_max)
				ax[i,j].set_ylim(-0.1, flat_max)
				if i<8:
					plt.setp( ax[i,j].get_xticklabels(), visible=False)
    			#if j>0:
    			plt.setp( ax[i,j].get_yticklabels(), visible=False)
		plt.tight_layout()
		plt.subplots_adjust(wspace=0, hspace=0)
		plt.savefig('../plots/flat/digplot.png', dpi=190)
		plt.show()
		'''
		plt.figure(1, (9,9))
		for i in range(9):
			for j in range(i):
				plt.plot(star_flats[i].flatten(), star_flats[j].flatten(), '.k')
				plt.xlabel('flat %d'%i)
				plt.ylabel('flat %d'%j)
				plt.xlim(-0.1, flat_max)
				plt.ylim(-0.1, flat_max)
				plt.tight_layout()
				plt.savefig('../plots/flat/flat_%d_%d.png'%(i,j), dpi=190)
				plt.clf()
		plt.plot(bkg_flats[0].flatten(), bkg_flats[1].flatten(), '.k')
		plt.show()
		m_flat = np.median(star_flats, axis=0)
		print m_flat.shape
		dev = np.square(star_flats-m_flat)
		m_dev = np.median(dev, axis=0)

		m_flat_bkg = np.median(bkg_flats, axis=0)
		dev_bkg = np.square(bkg_flats-m_flat_bkg)
		m_dev_bkg = np.mean(dev_bkg, axis=0)
		dev_median = np.median(m_dev_bkg)

		plt.imshow(m_dev, origin='lower', interpolation='none')
		plt.set_cmap('gray')
		plt.colorbar()
		plt.savefig('../plots/flat/var_map_normal100.png', dpi=190)
		#plt.show()

		plt.imshow(m_flat, origin='lower', interpolation='none')
		plt.set_cmap('gray')
		plt.colorbar()
		plt.savefig('../plots/flat/median_flat_normal100.png', dpi=190)
		#plt.show()

		plt.imshow(m_dev_bkg, origin='lower', interpolation='none', vmin=0.1*dev_median, vmax=10*dev_median)
		plt.set_cmap('gray')
		plt.colorbar()
		#plt.savefig('../plots/flat/var_map_bkg100.png', dpi=190)
		plt.show()

		plt.imshow(m_flat_bkg, origin='lower', interpolation='none')
		plt.set_cmap('gray')
		plt.colorbar()
		#plt.savefig('../plots/flat/median_bkg_flat100.png', dpi=190)
		plt.show()