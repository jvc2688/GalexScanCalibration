from astropy.io import fits as pyfits
import numpy as np
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import os
from astropy.convolution import convolve, Gaussian2DKernel, convolve_fft

if __name__ == '__main__':
    if False:
        print 'count'
        region='0-10'
        date='09-15-2017'
        hdu = pyfits.open('/scratch/dw1519/galex/fits/scan_map/'+date+'/count_map_'+region+'_count.fits')[0]
        img = hdu.data

        gauss = Gaussian2DKernel(stddev=5)
        im1 = convolve(img, gauss)

        hdu_w = pyfits.PrimaryHDU(im1)
        hdu_w.header = hdu.header
        hdulist = pyfits.HDUList([hdu_w])
        hdulist.writeto('/scratch/dw1519/galex/fits/scan_map/'+date+'/count_map_'+region+'_count_smooth.fits', clobber=False)

    if True:
        print 'exp'
        region='0-10'
        date='09-15-2017'
        hdu = pyfits.open('/scratch/dw1519/galex/fits/scan_map/'+date+'/count_map_'+region+'_exp.fits')[0]
        img = hdu.data

        gauss = Gaussian2DKernel(stddev=5)
        im1 = convolve(img, gauss)

        hdu_w = pyfits.PrimaryHDU(im1)
        hdu_w.header = hdu.header
        hdulist = pyfits.HDUList([hdu_w])
        hdulist.writeto('/scratch/dw1519/galex/fits/scan_map/'+date+'/count_map_'+region+'_exp_smooth.fits', clobber=False)