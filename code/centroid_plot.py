import c3
from astropy.io import fits as pyfits
import numpy as np
from astropy import wcs as pywcs
import matplotlib.pyplot as plt
import sys

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
		filename = '../fits/co/co_map%d_%d_zoom_large.fits'%(i,i+1)
		centroid = find_centroid(filename)
		centers.append(centroid)
	centers = np.array(centers)
	return centers

if __name__ == '__main__':
	initial = 300
	final = 1342
	centers = get_centers(initial, final)
	trange = np.arange(final-initial+1)+initial
	np.save('../data/offsets.npy', centers)
	'''
	plt.plot(trange, centers[:,0],'-b')
	plt.xlabel('time/s')
	plt.ylabel('RA/degree')
	plt.savefig('../plots/ra%d_%d_r.png'%(initial,final),dpi=190)
	plt.clf()
	plt.plot(trange, centers[:,1],'-b')
	plt.xlabel('time/s')
	plt.ylabel('DEC/degree')
	plt.savefig('../plots/dec%d_%d_r.png'%(initial,final),dpi=190)
	'''	
