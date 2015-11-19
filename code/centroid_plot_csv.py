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

if __name__ == '__main__':
	if True:
		scan_num = sys.argv[1]
		initial = int(sys.argv[2])#300
		final = int(sys.argv[3])#1343
		centers = []
		for i in range(initial, final+1):
			centers.append(np.load('../data/%s_new/cata/centroids%d.npy'%(scan_num, i)))
		centroids = np.concatenate(centers, axis=0)
		#centroids = -centroids
		print centroids.shape

		trange = np.arange((final-initial+1)*10)/10.+initial-0.5
		np.save('../data/%s_new/cata/offsets%d_10_new.npy'%(scan_num, initial), centroids)
		
		plt.plot(trange, centroids[:,0],'.b', markersize=2)
		plt.xlabel('time/s')
		plt.xlim(initial-2,final+2)
		#plt.ylim(-0.005, 0.005)
		plt.ylabel('RA/degree')
		plt.savefig('../plots/photon_%s_new/cata/ra%d_%d_new_r.png'%(scan_num,initial,final),dpi=190)
		plt.clf()
		plt.plot(trange, centroids[:,1],'.b', markersize=2)
		plt.xlabel('time/s')
		plt.xlim(initial-2,final+2)
		#plt.ylim(-0.005, 0.005)
		plt.ylabel('DEC/degree')
		plt.savefig('../plots/photon_%s_new/cata/dec%d_%d_new_r.png'%(scan_num,initial,final),dpi=190)
		plt.clf()
		
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
