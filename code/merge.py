import imagetools
import numpy as np
from astropy.io import fits as pyfits

if __name__ == '__main__':
	scan_list = ['05','14','23','32','41','50','59','68']

	tranges = [[0,0]]
 	skypos_f, skyrange_f = ([0, 0], [1.33333336,1.33333336])
	imsz_f = imagetools.deg2pix(skypos_f, skyrange_f, 0.0016666667)
 	count_f = np.zeros(imsz_f)
 	wcs_f = imagetools.define_wcs(skypos_f,skyrange_f,width=False,height=False,verbose=0,pixsz=0.0016666667)

	for scan_num in scan_list:
		hdu_list = pyfits.open('../fits/count_map_%s_gPr_cata_10_device_r0_nostar.fits'%scan_num)
		data = hdu_list[0].data
		#tranges.append(hdulist[0].header)
		count_f += data
	print np.median(count_f)*0.85, np.median(count_f)*1.15

	hdu = pyfits.PrimaryHDU(count_f)
	hdu = imagetools.fits_header('NUV', skypos_f, tranges, skyrange_f, hdu=hdu, wcs=wcs_f)
	hdulist = pyfits.HDUList([hdu])
	hdulist.writeto('../fits/count_map_05-68_gPr_cata_10_device_r0_nostar.fits', clobber=False)

