import numpy as np 
from astropy.io import fits as pyfits
import re
import glob

if __name__ == '__main__':

	if True:
		name_list = glob.glob("../name_new/p2/*")
		name_list = ['name2057-2102']
		name_list = ['name_185', 'name_239', 'name_248', 'name_221','name_230', 'name_3029', 'name_3038', 'name_3047', 'name_3245', 'name_3299', 'name_3317', 'name_3488']
		name_list = ['name_29s75', 'name_2984', 'name_3020']
		for name in name_list:
			#name = re.split('/', name)[3]
			if name == 'name2444-2480':
				continue
			print name
 			count = pyfits.open('../fits/scan_map/count_map_%s_gal_sec_count.fits'%name)
			exp = pyfits.open('../fits/scan_map/count_map_%s_gal_sec_exp.fits'%name)

			cm = count[0].data
			em = exp[0].data

			em[np.where(em<0.00001)]=1.
			print np.min(em)
			print np.max(em)

			ratio = cm/em
			print np.min(cm)
			print np.max(cm)

			print np.min(ratio)
			print np.max(ratio)

			hdu = pyfits.PrimaryHDU(ratio)
			hdu.header = count[0].header
			hdulist = pyfits.HDUList([hdu])
			hdulist.writeto('../fits/scan_map/count_map_%s_gal_sec_in.fits'%name, clobber=False)

	if False:
		count = pyfits.open('../fits/count_map_name671_gPr_cata_10_gal.fits')
		exp = pyfits.open('../fits/exp/count_map_name671_gPr_flip_gal_dead_sec.fits')

		cm = count[0].data
		em = exp[0].data

		em[np.where(em<0.00001)]=1.
		print np.min(em)
		print np.max(em)

		ratio = cm/em
		print np.min(cm)
		print np.max(cm)

		print np.min(ratio)
		print np.max(ratio)

		hdu = pyfits.PrimaryHDU(ratio)
		hdu.header = count[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/count_map_name671_gPr_cata_10_corr_gal.fits', clobber=False)

	#focal plane
	if False:
		focal = pyfits.open('count_map_05-68_gPr_cata_10_device_r0_nostar.fits')
		flat = pyfits.open('NUV_flat.fits')

		focal_m = focal[0].data
		flat_m = flat[0].data

		#flat_m = np.flipud(np.rot90(flat_m))
		
		flat_m = np.swapaxes(flat_m,0,1)

		print np.median(focal_m)*0.85, np.median(focal_m)*1.15
		print np.median(flat_m)*0.85, np.median(flat_m)*1.15

		#flat_m[np.where(flat_m<0.00001)]=1.
		flat_m[np.where(flat_m<0.001)]=1.

		#print np.min(flat_m)
		#print np.max(flat_m)

		ratio = focal_m/flat_m
		#print np.min(focal_m)
		#print np.max(focal_m)

		#print np.min(ratio)
		#print np.max(ratio)

		print np.median(ratio)*0.85, np.median(ratio)*1.15


		hdu = pyfits.PrimaryHDU(ratio)
		hdu.header = focal[0].header
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('count_map_05_68_gPr_cata_10_device_r0_ratio_hi.fits', clobber=False)

