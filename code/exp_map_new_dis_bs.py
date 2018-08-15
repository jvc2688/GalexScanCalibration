import numpy as np
from astropy.io import fits as pyfits
import imagetools
import gnomonic as gn
import pos_range
import scipy
import math
import threading
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from multiprocessing import Process
from multiprocessing import Manager
import sys
import os
import gc
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d



def exp_dead_new(file_num, name_file, imsz, wcs, hrflat, asp_solution, dead, return_dict):
	count = np.zeros(imsz)
	rotate = scipy.ndimage.interpolation.rotate

	x_lim = imsz[0]
	y_lim = imsz[1]

	length = hrflat.shape[0]
	half_len = length/2.
	print half_len
	coo = asp_solution[:,1:3]
	foc_list = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
	print asp_solution.shape[0]
	for i in range(asp_solution.shape[0]):
		res = rotate(hrflat,-asp_solution[i,3],reshape=False,order=1,prefilter=False)
		foc = foc_list[i,:]#wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
		cent = np.zeros(2)
		cent[0] = np.floor(foc[1]-0.5)
		cent[1] = np.floor(foc[0]-0.5)
		#print cent
		res_c_x = res.shape[0]/2.
		res_c_y = res.shape[1]/2.

		xl = int(cent[0]-half_len)
		xh = int(cent[0]+half_len)

		yl = int(cent[1]-half_len)
		yh = int(cent[1]+half_len)

		if xl<=0:
			dxl = -xl
			xl = 0
		else:
			dxl = 0

		if xh>x_lim:
			dxh = xh-x_lim
			xh = x_lim
		else:
			dxh = 0

		if yl<=0:
			dyl = -yl
			yl = 0
		else:
			dyl = 0

		if yh>y_lim:
			dyh = yh-y_lim
			yh = y_lim
		else:
			dyh = 0

		#res_x, res_y = res.shape[0], res.shape[1]
		res_xl = int(dxl)
		res_xh = int(length-dxh)
		res_yl = int(dyl)
		res_yh = int(length-dyh)
		count[xl:xh, yl:yh] += res[res_xl:res_xh, res_yl:res_yh]*0.02*(1-dead[i])

		#print cent[0], cent[1]
		#print xl,xh,yl,xh
		#count[xl:xh, yl:yh] += res*0.005*(1-dead[i])#res[res_c_x-half_len:res_c_x+half_len, res_c_y-half_len:res_c_y+half_len]*0.005*(1-dead[i])
		if i%2000==0:
			with open('../fits/scan_map/%s_gal_sec_exp_tmp%d.dat'%(name_file, file_num),'w') as f:
				f.write('%d'%i)
			print i
	print '%d done'%file_num
	#return_dict[file_num] = count
	np.save('../fits/scan_map/%s_gal_sec_exp_tmp%d.npy'%(name_file, file_num), count)



if __name__ == '__main__':

#process
	if True:
		name_file = sys.argv[1]
		with open('../name_scan/%s'%name_file) as f:
			name_list = f.read().splitlines()
		print name_list

		if len(sys.argv)>3:
			suffix = sys.argv[3]
		else:
			suffix = ''

		num_t = int(sys.argv[2])
		pixsz = 0.0005555555556 #0.000416666666666667 #0.00166667
		asp_solution_list = []
		ix_list = []
		dead_list = []
		corr_list = []
		gal_l = []
		for name in name_list:
			try:
				asp_solution_tmp = np.load('../data/photon_list/{0}{1}-new-fine_asp_new_bs.npy'.format(name, suffix))
			except IOError:
				co_data = pyfits.open('../AIS_GAL_SCAN/asprta/{0}{1}-new-fine-asprta.fits'.format(name, suffix))[1].data
				T = co_data['T']
				ra = co_data['ra']
				dec = co_data['dec']
				roll = co_data['roll']
				n = 5
				t_new = np.arange((T.shape[0]-1)*n)*0.02+T[0]

				'''
				f_ra = interp1d(T, ra, kind='cubic')
				f_dec = interp1d(T, dec, kind='cubic')
				f_roll = interp1d(T, roll, kind='cubic')

				ra_new = f_ra(t_new)
				dec_new = f_dec(t_new)
				roll_new = f_roll(t_new)
				'''
				
				spl_ra = splrep(T, ra)
				spl_dec = splrep(T, dec)
				spl_roll = splrep(T, roll)
				ra_new = splev(t_new, spl_ra)
				dec_new = splev(t_new, spl_dec)
				roll_new = splev(t_new, spl_roll)
			
				'''
				ra_new = np.interp(t_new, T, ra)
				dec_new = np.interp(t_new, T, dec)
				roll_new = np.interp(t_new, T, roll)
				'''
				asp_solution_tmp = np.array([t_new, ra_new, dec_new, roll_new]).T
				np.save('../data/photon_list/{0}{1}-new-fine_asp_new_bs.npy'.format(name, suffix), asp_solution_tmp)
			#ix_tmp = (np.round(asp_solution_tmp[:,0]-asp_solution_tmp[0,0])+1).astype(int)
			scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%name)[1].data
			scst_time = scst_tmp['pktime']
			ix_tmp = np.digitize(asp_solution_tmp[:,0], scst_time) 
			ix_mask = ix_tmp<scst_time.shape[0]
			ix_tmp = ix_tmp[ix_mask]
			asp_solution_tmp = asp_solution_tmp[ix_mask]
			hv = scst_tmp['hvnom_nuv']
			mask = hv[ix_tmp]>0
			ix_tmp = ix_tmp[mask]
			asp_solution_tmp = asp_solution_tmp[mask]
			asp_solution_list.append(asp_solution_tmp)			
			#limit = scst_tmp['t_dead_nuv'].shape[0]-1
			#ix_tmp[ix_tmp>limit] = limit
			dead = scst_tmp['t_dead_nuv']
			tec = scst_tmp['NDCTEC']
			fec = scst_tmp['NDCFEC']
			fec_mask = fec>150000
			width=1000
			while True:
				ratio_mask = (fec>(150000-width)) & (fec<(150000+width))
				if np.sum(ratio_mask)>10:
					break
				else:
					width+=1000
			ratio = np.mean((1.-dead[ratio_mask])/(tec[ratio_mask]/fec[ratio_mask]))
			print 'ratio:{0}'.format(ratio)
			dead[fec_mask] = 1.-ratio*tec[fec_mask]/fec[fec_mask]
			dead_tmp = dead[ix_tmp]
			dead_list.append(dead_tmp)

			#asp_old = np.load('../data/photon_list/%s_asp.npy'%name)
			#sky_data = SkyCoord(asp_old[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
			#calculate the image size and location
			asp_map = np.load('../data/photon_list/{0}{1}_asp.npy'.format(name, suffix))
			sky_data = SkyCoord(asp_map[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
			gal = sky_data.transform_to(Galactic)
			asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
			gal_l.append(np.mean(asprta[:,0]))

		asp_solution = np.concatenate(asp_solution_list, axis=0)
		sky_data = SkyCoord(asp_solution[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
		gal = sky_data.transform_to(Galactic)
		asp_solution[:,1:3] = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
	
		dead = np.concatenate(dead_list, axis=0)

		asp_solutions = np.array_split(asp_solution, num_t)
		deads = np.array_split(dead, num_t)

		flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
		flat =flat_hdu[0].data
		size=flat_hdu[0].header['CDELT2']
		resample=3 #4
		hrflat = scipy.ndimage.interpolation.zoom(flat, resample, order=0,prefilter=False)#scipy.ndimage.interpolation.zoom(flat,size/pixsz,order=1,prefilter=False)
		print hrflat.shape
		hrflat = np.swapaxes(hrflat,0,1)
		print hrflat.shape
		hrflat = scipy.ndimage.interpolation.rotate(hrflat,-23.4,reshape=False,order=1,prefilter=False)
		print hrflat.shape
		skypos = [np.mean(gal_l), 0.]
		skyrange = [np.max(gal_l)-np.min(gal_l)+1.5, 21.5]
		tranges = []
		trange = [0, 0]
		imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz).astype(int)
		wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)

		manager = Manager()
		return_dict = manager.dict()
		p_list = []

		for file_num in range(0, num_t):
			p = Process(target=exp_dead_new, args=(file_num, name_file, imsz, wcs, hrflat, asp_solutions[file_num], deads[file_num], return_dict))
			p.start()
			p_list.append(p)

		for p in p_list:
			p.join()
		print 'all done'
		
		count = np.zeros(imsz)

		for i in range(0, num_t):
			count += np.load('../fits/scan_map/%s_gal_sec_exp_tmp%d.npy'%(name_file, i)) #return_dict[i]
			os.remove('../fits/scan_map/%s_gal_sec_exp_tmp%d.npy'%(name_file, i))

		tranges.append(trange)
		hdu = pyfits.PrimaryHDU(count)
		hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/scan_map/count_map_{0}{1}_exp_10.fits'.format(name_file, suffix), clobber=False)


#process, small pitach
	if False:
		name_file = sys.argv[1]
		with open('../name_new/%s'%name_file) as f:
			name_list = f.read().splitlines()
		print name_list

		num_t = int(sys.argv[2])
		asp_solution_list = []
		ix_list = []
		dead_list = []
		gal_l = []
		for name in name_list:
			asp_solution_tmp = np.load('../data/photon_list/%s_asp_cal.npy'%name)
			asp_solution_list.append(asp_solution_tmp)			
			ix_tmp = (np.round(asp_solution_tmp[:,0]-asp_solution_tmp[0,0])+1).astype(int)
			dead_tmp = np.zeros(ix_tmp.shape[0])
			scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%name)[1].data
			limit = scst_tmp['t_dead_nuv'].shape[0]-1
			ix_tmp[ix_tmp>limit] = limit
			dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]
			dead_list.append(dead_tmp)

			#asp_old = np.load('../data/photon_list/%s_asp.npy'%name)
			#sky_data = SkyCoord(asp_old[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
			sky_data = SkyCoord(asp_solution_tmp[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
			gal = sky_data.transform_to(Galactic)
			asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
			gal_l.append(np.mean(asprta[:,0]))

		asp_solution = np.concatenate(asp_solution_list, axis=0)

		sky_data = SkyCoord(asp_solution[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
		gal = sky_data.transform_to(Galactic)
		asp_solution[:,1:3] = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
	
		dead = np.concatenate(dead_list, axis=0)

		asp_solutions = np.array_split(asp_solution, num_t)
		deads = np.array_split(dead, num_t)
		flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
		flat =flat_hdu[0].data
		size=flat_hdu[0].header['CDELT2']
		hrflat = scipy.ndimage.interpolation.zoom(flat,size/0.002,order=1,prefilter=False)
		print hrflat.shape
		hrflat = np.swapaxes(hrflat,0,1)
		print hrflat.shape
		hrflat = scipy.ndimage.interpolation.rotate(hrflat,59.,reshape=False,order=1,prefilter=False)
		print hrflat.shape
		skypos = [np.mean(gal_l), 0.]
		skyrange = [np.max(gal_l)-np.min(gal_l)+2., 24.]
		tranges = []
		trange = [0, 0]
		imsz = imagetools.deg2pix_g(skypos, skyrange, 0.002)
		wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)

		manager = Manager()
		return_dict = manager.dict()
		p_list = []

		for file_num in range(0, num_t):
			p = Process(target=exp_dead_new, args=(file_num, imsz, wcs, hrflat, asp_solutions[file_num], deads[file_num], return_dict))
			p.start()
			p_list.append(p)

		for p in p_list:
			p.join()
		print 'all done'
		
		count = np.zeros(imsz)

		for i in range(0, num_t):
			count += return_dict[i]

		tranges.append(trange)
		hdu = pyfits.PrimaryHDU(count)
		hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/count_map_%s_gPr_flip_gal_dead.fits'%name_file, clobber=False)


