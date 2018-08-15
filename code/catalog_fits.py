from __future__ import print_function
import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc
from astropy.io import fits as pyfits
import astropy.coordinates as pycoo

def get_catalog_t(skypos, angle):
	center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
	rad = np.cos(angle*np.pi/180.)
	#print(center, rad)
	#load catalog
	hdulist = pyfits.open('../data/bstar.fits')
	length = int(hdulist[1].header['NAXIS2'])
	data = hdulist[1].data
	stars = np.zeros((length,2))
	stars[:,0] = data['ra']
	stars[:,1] = data['dec']
	stars_rad = stars*np.pi/180.
	#convert (ra, dec) to (x, y, z)
	X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
	stars_car = np.array([X,Y,Z], dtype='float64').T

	sep = np.dot(stars_car, center)
	coo = stars[sep>=rad, :]

	coo_list = []
	for i in range(0,coo.shape[0]):
		tmp_coo = (coo[i,0], coo[i,1])
		coo_list.append(tmp_coo)

	return coo_list


def get_catalog(skypos, angle):
	center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
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
	X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
	stars_car = np.array([X,Y,Z], dtype='float64').T

	sep = np.dot(stars_car, center)
	coo = stars[sep>=rad, :]

	return coo

def get_catalog_tgas(skypos, angle):
	center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
	rad = np.cos(angle*np.pi/180.)
	print(center, rad)
	#load catalog
	hdulist = pyfits.open('../data/bstar.fits')
	length = int(hdulist[1].header['NAXIS2'])
	data = hdulist[1].data
	data = np.load('../data/tgas.npy')
	length = data.shape[0]
	stars = np.zeros((length,2))
	stars[:,0] = data['ra']
	stars[:,1] = data['dec']
	stars_rad = stars*np.pi/180.
	#convert (ra, dec) to (x, y, z)
	X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
	stars_car = np.array([X,Y,Z], dtype='float64').T

	sep = np.dot(stars_car, center)
	coo = stars[sep>=rad, :]

	return coo

def get_catalog_tycho(skypos, angle):
	center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
	rad = np.cos(angle*np.pi/180.)
	print(center, rad)
	#load catalog
	hdulist = pyfits.open('../data/tycho2.fits')
	length = int(hdulist[1].header['NAXIS2'])
	data = hdulist[1].data
	mask = data['BTmag']>8
	data = data[mask]
	length = data.shape[0]
	stars = np.zeros((length,2))
	stars[:,0] = data['RAJ2000']
	stars[:,1] = data['DEJ2000']
	stars_rad = stars*np.pi/180.
	#convert (ra, dec) to (x, y, z)
	X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
	stars_car = np.array([X,Y,Z], dtype='float64').T

	sep = np.dot(stars_car, center)
	coo = stars[sep>=rad, :]

	return coo

def get_catalog_tycho_fast(skypos, angle, stars, stars_car, data):
	center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
	rad = np.cos(angle*np.pi/180.)
	print(center, rad)
	mask = (data['BTmag']>8) #& (data['BTmag']<18)
	sep = np.dot(stars_car, center)
	mask2 = sep>=rad
	coo = stars[mask&mask2, :]

	return coo

def get_catalog_tycho_square(up,down,left,right):
	hdulist = pyfits.open('../data/tycho2.fits')
	length = int(hdulist[1].header['NAXIS2'])
	data = hdulist[1].data
	mask = data['BTmag']>8
	data = data[mask]
	length = data.shape[0]
	stars = np.zeros((length,2))
	stars[:,0] = data['Glon']
	stars[:,1] = data['Glat']
	mask = (stars[:,0]<=left) & (stars[:,0]>=right) & (stars[:,1]<=down) & (stars[:,1]>=up)
	#stars_rad = stars*np.pi/180.
	#convert (ra, dec) to (x, y, z)
	#X, Y, Z = pycoo.spherical_to_cartesian(1, stars_rad[:,1], stars_rad[:,0])
	#stars_car = np.array([X,Y,Z], dtype='float64').T

	#sep = np.dot(stars_car, center)
	coo = stars[mask, :]

	return coo

def get_pos_time(filename, time):
	hdulist = pyfits.open(filename)
	co_data = hdulist[1].data
	intitial_asp = co_data[time]
	skypos = np.array([intitial_asp[1], intitial_asp[2]])
	intitial_time = intitial_asp[0]
	return skypos, intitial_time

def get_catalog_matched(skypos, angle, stars, stars_car, data):
	center = pycoo.spherical_to_cartesian(1, skypos[1]*np.pi/180., skypos[0]*np.pi/180.)
	rad = np.cos(angle*np.pi/180.)
	#print(center, rad)
	#mask = (data['BTmag']>8) #& (data['BTmag']<18)
	sep = np.dot(stars_car, center)
	mask = sep>=rad
	coo = stars[mask, :]

	return coo, data[mask]

if __name__ == '__main__':
  intitial_sec = 451
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
  hdulist.writeto('../fits/photons_cata/catalog_5_%d_nom_try.fits'%(intitial_sec), clobber=False)
