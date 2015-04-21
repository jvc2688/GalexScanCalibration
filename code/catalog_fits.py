from __future__ import print_function
import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc
from astropy.io import fits as pyfits
import astropy.coordinates.distances as dis

def get_catalog(skypos, angle):
	center = dis.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
	rad = np.cos(angle*np.pi/180.)
	print(center, rad)
	#load catalog
	hdulist = pyfits.open('../data/bstar.fits')
	length = int(hdulist[1].header['NAXIS2'])
	data = hdulist[1].data
	stars = np.zeros((length,2))
	stars[:,0] = data['ra']
	stars[:,1] = data['dec']
	stars_rad = stars*np.pi/180.
	#convert (ra, dec) to (x, y, z)
	X, Y, Z = dis.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
  	stars_car = np.array([X,Y,Z], dtype='float64').T

  	sep = np.dot(stars_car, center)
  	coo = stars[sep>=rad, :]

  	return coo

def get_pos_time(filename, time):
	hdulist = pyfits.open(filename)
  	co_data = hdulist[1].data
  	intitial_asp = co_data[time]
  	skypos = np.array([intitial_asp[1], intitial_asp[2]])
  	intitial_time = intitial_asp[0]
  	return skypos, intitial_time

if __name__ == '__main__':
	intitial_sec = 100
  	skypos, intitial_time = get_pos_time('../data/AIS_GAL_SCAN_00005_0001-asprta.fits', intitial_sec)
  	coo = get_catalog(skypos, 0.69)
  	skyrange = [1.38, 1.38]

  	print(coo.shape)

  	imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
  	count = np.zeros(imsz)

  	wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
	foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
	H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                                   bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
	count += H

	tranges = [[0,1000]]
  	hdu = pyfits.PrimaryHDU(count)
  	hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu)
  	hdulist = pyfits.HDUList([hdu])
  	hdulist.writeto('catalog_5_%d_nom_try.fits'%(intitial_sec), clobber=False)
