import c3
from astropy.io import fits as pyfits
import numpy as np
from astropy import wcs as pywcs
import matplotlib.pyplot as plt
import sys
import aplpy
import os

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

def corr_plot(filename, title):
	print filename
	centroid = find_centroid(filename)
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

if __name__ == '__main__':
	if False:
		initial = 300
		final = 1342
		centers = get_centers(initial, final)
		trange = np.arange(final-initial+1)+initial
		np.save('../data/offsets300-1342_r.npy', centers)
		
		plt.plot(trange, centers[:,0],'-b')
		plt.xlabel('time/s')
		plt.ylabel('RA/degree')
		plt.savefig('../plots/ra%d_%d_new.png'%(initial,final),dpi=190)
		plt.clf()
		plt.plot(trange, centers[:,1],'-b')
		plt.xlabel('time/s')
		plt.ylabel('DEC/degree')
		plt.savefig('../plots/dec%d_%d_new.png'%(initial,final),dpi=190)

	if False:
		initial = 600
		final = 1342
		for i in range(initial, final+1):
			filename = '../fits/co/right/co_map%d_%d_zoom_nocal.fits'%(i,i+1)
			corr_plot(filename, 'corr between %ds and %ds'%(i,i+1))

	if False:
		initial = 500
		final = 510
		centers = get_centers_half(initial, final)

		trange = np.arange((final-initial+1)*10)/10.+initial-0.5
		np.save('../data/offsets300-899_half_r.npy', centers)
		
		plt.plot(trange, centers[:,0],'.b')
		plt.xlabel('time/s')
		plt.ylabel('RA/degree')
		plt.savefig('../plots/ra%d_%d_half_new.png'%(initial,final),dpi=190)
		plt.clf()
		plt.plot(trange, centers[:,1],'.b')
		plt.xlabel('time/s')
		plt.ylabel('DEC/degree')
		plt.savefig('../plots/dec%d_%d_half_new.png'%(initial,final),dpi=190)
		
	if False:
		initial = 500
		final = 510
		for i in range(initial, final+1):
			for j in range(10):
				filename = '../fits/test/co_map%d_%d_10.fits'%(i,j)
				corr_plot(filename, 'corr between %.1fs and %.1fs'%(i+0.1*j-0.5,i+0.1*(j+1)-0.5))

	if True:
		initial = 500
		final = 510

		trange = np.arange((final-initial+1)*10)/10.+initial-0.5
		centers = np.load('../data/offsets300-899_half_r.npy')

		for i in range(1, 110):
			centers[i] = centers[i]+centers[i-1]
		
		plt.clf()
		plt.plot(trange, centers[:,0],'.b')
		plt.xlabel('time/s')
		plt.ylabel('RA/degree')
		plt.savefig('../plots/ra%d_%d_half_cul.png'%(initial,final),dpi=190)
		plt.clf()
		plt.plot(trange, centers[:,1],'.b')
		plt.xlabel('time/s')
		plt.ylabel('DEC/degree')
		plt.savefig('../plots/dec%d_%d_half_cul.png'%(initial,final),dpi=190)
