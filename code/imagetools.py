"""
This file is part of the Galex Photon Catalog project.
Copyright 2012 Chase Million
"""

import numpy as np
from astropy import wcs as pywcs
from astropy.io import fits as pyfits
import scipy.misc
import scipy.special # erfc
import scipy.ndimage

def deg2pix(skypos,skyrange,pixsz=0.000416666666666667):
	"""Converts degrees to GALEX pixels rounded up to the nearest pixel
	so that the number of degrees specified will fully fit into the frame.
	"""
    
	wcs = pywcs.WCS(naxis=2) # NAXIS = 2
	wcs.wcs.cdelt = [pixsz,pixsz]
	wcs.wcs.ctype = ['RA---TAN','DEC--TAN']
	# Set the reference pixel to [0,0] (which is actually the default)
	wcs.wcs.crpix = [0.,0.]
	# Fix the lower bound [RA, Dec] to the reference pixel
	wcs.wcs.crval = [skypos[0]-skyrange[0]/2.,skypos[1]-skyrange[1]/2.]
	# Find the pixel position of the upper bound [RA,Dec]
	coo = [skypos[0]+skyrange[0]/2.,skypos[1]+skyrange[1]/2.]
	
	return np.abs(np.floor(wcs.sip_pix2foc(wcs.wcs_world2pix([coo],1),1)[0]))[::-1]

def define_wcs(skypos,skyrange,width=False,height=False,verbose=0,
               pixsz=0.002):
	"""Define the world coordinate system (WCS)."""
	wcs = pywcs.WCS(naxis=2) # NAXIS = 2
	imsz = deg2pix(skypos,skyrange, pixsz)#, 0.002)
	wcs.wcs.cdelt = np.array([-pixsz,pixsz])#wcs.wcs.cdelt = np.array([-pixsz,pixsz])
	wcs.wcs.ctype = ['RA---TAN','DEC--TAN']
	wcs.wcs.crpix = [(imsz[1]/2.)+0.5,(imsz[0]/2.)+0.5]
	wcs.wcs.crval = skypos
	return wcs

def deg2pix_g(skypos,skyrange,pixsz=0.000416666666666667):
	"""Converts degrees to GALEX pixels rounded up to the nearest pixel
	so that the number of degrees specified will fully fit into the frame.
	"""
    
	wcs = pywcs.WCS(naxis=2) # NAXIS = 2
	wcs.wcs.cdelt = [pixsz,pixsz]
	wcs.wcs.ctype = ['GLON-TAN','GLAT-TAN']
	# Set the reference pixel to [0,0] (which is actually the default)
	wcs.wcs.crpix = [0.,0.]
	# Fix the lower bound [RA, Dec] to the reference pixel
	wcs.wcs.crval = [skypos[0]-skyrange[0]/2.,skypos[1]-skyrange[1]/2.]
	# Find the pixel position of the upper bound [RA,Dec]
	coo = [skypos[0]+skyrange[0]/2.,skypos[1]+skyrange[1]/2.]
	
	return np.abs(np.floor(wcs.sip_pix2foc(wcs.wcs_world2pix([coo],1),1)[0]))[::-1]

def define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,
               pixsz=0.002):
	"""Define the world coordinate system (WCS)."""
	wcs = pywcs.WCS(naxis=2) # NAXIS = 2
	imsz = deg2pix_g(skypos,skyrange, pixsz)#, 0.002)
	wcs.wcs.cdelt = np.array([-pixsz,pixsz])#wcs.wcs.cdelt = np.array([-pixsz,pixsz])
	wcs.wcs.ctype = ['GLON-TAN','GLAT-TAN']
	wcs.wcs.crpix = [(imsz[1]/2.)+0.5,(imsz[0]/2.)+0.5]
	wcs.wcs.crval = skypos
	return wcs

def fits_header(band,skypos,tranges,skyrange,width=False,height=False,
				verbose=0,tscale=1000.,hdu=False,wcs=None,retries=20):
	"""Populate a FITS header."""
	hdu = hdu if hdu else pyfits.PrimaryHDU()
	if wcs==None:
		wcs = define_wcs(skypos,skyrange,width=width,height=height)
	hdu.header['CDELT1'],hdu.header['CDELT2'] = wcs.wcs.cdelt
	hdu.header['CTYPE1'],hdu.header['CTYPE2'] = wcs.wcs.ctype
	hdu.header['CRPIX1'],hdu.header['CRPIX2'] = wcs.wcs.crpix
	hdu.header['CRVAL1'],hdu.header['CRVAL2'] = wcs.wcs.crval
	#hdu.header['RA_CENT'],hdu.header['DEC_CENT'] = wcs.wcs.crval # Dupe.
	hdu.header['EQUINOX'],hdu.header['EPOCH'] = 2000., 2000.
	hdu.header['BAND'] = 1 if band=='NUV' else 2

	# Put the total exposure time into the primary header
	hdu.header['EXPTIME'] = 0.
	for trange in tranges:
		hdu.header['EXPTIME'] += trange[1]-trange[0] #dbt.compute_exptime(band,trange,verbose=verbose,retries=retries)

	if len(tranges)==1:
	# Put the time range into the primary header for a single frame image
		hdu.header['EXPSTART'],hdu.header['EXPEND'] = tranges[0]
		# These are the proper keywords for this:
		hdu.header['TIME-OBS'],hdu.header['TIME-END'] = tranges[0]

	return hdu


