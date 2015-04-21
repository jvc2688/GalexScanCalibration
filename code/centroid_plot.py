import c3
from astropy.io import fits as pyfits
import numpy as np
from astropy import wcs as pywcs
import matplotlib.pyplot as plt

if __name__ == '__main__':
	centers = []
	for i in range(1,91):
		hdulist = pyfits.open('../fits/co/co_map%d_%d_zoom_large.fits'%(i,i+1))
		w = pywcs.WCS(hdulist[0].header, hdulist)
		data = hdulist[0].data
		data = data.byteswap(True).newbyteorder()
		cy, cx = c3.find_centroid(data)
		centroid = w.wcs_pix2world(w.sip_foc2pix([[cx, cy]],1),1)[0]
		if centroid[0]>1:
			centroid[0] = centroid[0]-360.
		centers.append(centroid)
	centers = np.array(centers)
	plt.plot(centers[:,0],'.b')
	plt.xlabel('time/s')
	plt.ylabel('RA/degree')
	plt.savefig('../plots/ra1_90.png',dpi=190)
	plt.clf()
	plt.plot(centers[:,1],'.b')
	plt.xlabel('time/s')
	plt.ylabel('DEC/degree')
	plt.savefig('../plots/dec1_90.png',dpi=190)
