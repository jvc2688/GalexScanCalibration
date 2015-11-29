import numpy as np
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import pylab as P
import spilt_csv as spi
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
import imagetools
import matplotlib.pyplot as plt

def get_star_mask(scan_name, fits=False):
	num_t = 1
	asp_solution_list = []
	ix_list = []
	dead_list = []
	gal_l = []
	
	asp_solution = np.load('../data/photon_list/%s_asp_cal.npy'%name)
	sky_data = SkyCoord(asp_solution[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
	gal = sky_data.transform_to(Galactic)
	asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
	gal_l.append(np.mean(asprta[:,0]))

	sky_data = SkyCoord(asp_solution[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
	gal = sky_data.transform_to(Galactic)
	asp_solution[:,1:3] = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)

	asp_solutions = np.array_split(asp_solution, num_t)
	skypos = [np.mean(gal_l), 0.]
	skyrange = [np.max(gal_l)-np.min(gal_l)+2., 24.]
	tranges = []
	trange = [0, 0]
	imsz = imagetools.deg2pix_g(skypos, skyrange, 0.00166667)
	wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.00166667)

	count = np.zeros(imsz)

	cata = spi.load_obj('../data/%s_starset_extra_all_star'%scan_name)
	cata_a = np.array(list(cata))
	cata_len = len(cata_a)
	print cata_a

	mask_faint = cata_a[:,3]<40
	mask_mid = np.logical_and(cata_a[:,3]>=40, cata_a[:,3]<120)
	mask_bright = np.logical_and(cata_a[:,3]>=120, cata_a[:,3]<250)
	mask_ultra = cata_a[:,3]>=250
	mask_list = [mask_faint, mask_mid, mask_bright, mask_ultra]
	size_list = [[-1,0,1],[-3,-2,-1,0,1,2,3],[-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9], 
							[-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]]

	for bright in range(4):
		star_mask = mask_list[bright]
		size = size_list[bright]
		cata_tmp = cata_a[star_mask]
		sky_data = SkyCoord(cata_tmp[:,0:2], unit='deg', frame=FK5, equinox='J2000.0')
		gal = sky_data.transform_to(Galactic)
		star_co = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)

		foc_list = wcs.sip_pix2foc(wcs.wcs_world2pix(star_co,1),1)-0.5
		print foc_list
		print np.max(foc_list[:,0]), np.max(foc_list[:,1])
		for i in size:
			for j in size:
				coo_tmp = np.zeros(foc_list.shape)
				coo_tmp[:,0] = foc_list[:,0]+i
				coo_tmp[:,1] = foc_list[:,1]+j
				mask_0 = np.logical_and(coo_tmp[:,0]>=0, coo_tmp[:,0]<1285)
				mask_1 = np.logical_and(coo_tmp[:,1]>=0, coo_tmp[:,1]<15310)
				mask = np.logical_and(mask_0, mask_1)
				coo_tmp = coo_tmp[mask]
				coo_tmp = np.array(coo_tmp, dtype=int)
				count[coo_tmp[:,1], coo_tmp[:,0]] += 1

	if fits:
		tranges = [[0,0]]
		hdu = pyfits.PrimaryHDU(count)
		hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/count_map_%s_gPr_flip_gal_star.fits'%scan_name, clobber=False)

	return count



if __name__ == '__main__':
	if True:

		scan_name = 'AIS_GAL_SCAN_00014_0001'

		star_mask = star_mask(scan_name, False)

		asp = np.load('../data/photon_list/%s_asp_cal.npy'%scan_name)
		asp_uni, uni_index=np.unique(asp[:,0], True)
		asp_uni = asp[uni_index]

		sky_data = SkyCoord(asp_uni[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
		gal = sky_data.transform_to(Galactic)
		asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)

		print asprta.shape

		in_map = pyfits.open('../fits/count_map_name14_gal_in.fits')[0].data
		print in_map.shape
		time_list = []
		bkg_list = []
		for i in range(2000, asprta.shape[0]/2-2000, 1000):
			time = asp_uni[i,0]
			print time
			start = wcs.sip_pix2foc(wcs.wcs_world2pix([asprta[i-1000]],1),1)-0.5
			end = wcs.sip_pix2foc(wcs.wcs_world2pix([asprta[i+1000]],1),1)-0.5
			in_region = in_map[end[0,1]:start[0,1], 300:1000].flatten()
			mask_region = star_mask[end[0,1]:start[0,1], 300:1000].flatten()
			bkg = np.mean(in_region[mask_region<1])
			#print bkg
			time_list.append(time)
			bkg_list.append(bkg)

		plt.plot(time_list, bkg_list, '.k')
		plt.xlabel('Time s')
		plt.ylabel('BKG photons/pixel/s')
		plt.savefig('../plots/bkg.png', dpi=190)



	#hist
	if False:
		size = 200
		data = np.zeros([size,size])
		star = np.zeros([size,size])
		for i in range(80):
			data += pyfits.open('../fits/flat/chunk/downsample_new/data_chunk%d.fits'%i)[0].data
			star += pyfits.open('../fits/flat/chunk/downsample_new/star_count%d.fits'%i)[0].data

		data = data.flatten()
		star = star.flatten()

		P.figure()
		n, bins, patches = P.hist([data,star], 70, histtype='bar',label=['data', 'star'])
		P.legend()
		plt.ylim(0,4000)
		plt.xlabel('number of photons')
		plt.ylabel('number of pixels')
		plt.show()
