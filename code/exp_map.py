import numpy as np
from astropy.io import fits as pyfits
import imagetools
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

class exp(threading.Thread):
    def __init__(self, thread_id, imsz, wcs, hrflat, asp_solution):
        threading.Thread.__init__(self)
        self.thread_id = thread_id
        self.count = np.zeros(imsz)
        self.wcs = wcs
        self.hrflat = hrflat
        self.asp_solution = asp_solution
    def run(self):
		print '%d start'%self.thread_id
		for i in range(self.asp_solution.shape[0]):
			res = scipy.ndimage.interpolation.rotate(self.hrflat,-self.asp_solution[i,3],reshape=True,order=1,prefilter=False)
			coo = [self.asp_solution[i,1:3]]
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
			cent = np.zeros(2)
			cent[0] = np.floor(foc[0]-0.5)[1]
			cent[1] = np.floor(foc[0]-0.5)[0]
			#print cent
			xl = cent[0]-res.shape[0]/2.
			dxl = xl-0 if xl-0<=0 else 0
			xl = xl if xl>=0 else 0
			xh = cent[0]+res.shape[0]/2.
			dxh = imsz[0]-xh if imsz[0]-xh<=0 else 0
			xh = xh if xh<imsz[0] else imsz[0]

			yl = cent[1]-res.shape[1]/2.
			dyl = yl-0 if yl-0<=0 else 0
			yl = yl if yl>=0 else 0
			yh = cent[1]+res.shape[1]/2.
			dyh = imsz[1]-yh if imsz[1]-yh<=0 else 0
			yh = yh if yh<imsz[1] else imsz[1]
			#print xl,xh,yl,yh, res.shape
			res_x, res_y = res.shape[0], res.shape[1]
			self.count[xl:xh, yl:yh] += res[-dxl:np.ceil(res_x+dxh), -dyl:np.ceil(res_y+dyh)]*0.005
			if i%2000==0:
				print i
		print '%d done'%self.thread_id

class exp_dead(threading.Thread):
    def __init__(self, thread_id, imsz, wcs, hrflat, asp_solution, dead):
        threading.Thread.__init__(self)
        self.thread_id = thread_id
        self.count = np.zeros(imsz)
        self.wcs = wcs
        self.hrflat = hrflat
        self.asp_solution = asp_solution
        self.dead = dead
    def run(self):
		print '%d start'%self.thread_id, self.dead.shape[0], self.asp_solution.shape[0]
		for i in range(self.asp_solution.shape[0]):
			res = scipy.ndimage.interpolation.rotate(self.hrflat,-self.asp_solution[i,3],reshape=True,order=1,prefilter=False)
			coo = [self.asp_solution[i,1:3]]
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
			cent = np.zeros(2)
			cent[0] = np.floor(foc[0]-0.5)[1]
			cent[1] = np.floor(foc[0]-0.5)[0]
			#print cent
			xl = cent[0]-res.shape[0]/2.
			dxl = xl-0 if xl-0<=0 else 0
			xl = xl if xl>=0 else 0
			xh = cent[0]+res.shape[0]/2.
			dxh = imsz[0]-xh if imsz[0]-xh<=0 else 0
			xh = xh if xh<imsz[0] else imsz[0]

			yl = cent[1]-res.shape[1]/2.
			dyl = yl-0 if yl-0<=0 else 0
			yl = yl if yl>=0 else 0
			yh = cent[1]+res.shape[1]/2.
			dyh = imsz[1]-yh if imsz[1]-yh<=0 else 0
			yh = yh if yh<imsz[1] else imsz[1]
			#print xl,xh,yl,yh, res.shape
			res_x, res_y = res.shape[0], res.shape[1]
			self.count[xl:xh, yl:yh] += res[-dxl:np.ceil(res_x+dxh), -dyl:np.ceil(res_y+dyh)]*0.005*(1-self.dead[i])
			if i%2000==0:
				print i
		print '%d done'%self.thread_id

class exp_dead_g(threading.Thread):
    def __init__(self, thread_id, imsz, wcs, hrflat, asp_solution, dead):
        threading.Thread.__init__(self)
        self.thread_id = thread_id
        self.count = np.zeros(imsz)
        self.wcs = wcs
        self.hrflat = hrflat
        self.asp_solution = asp_solution
        self.dead = dead
    def run(self):
		print '%d start'%self.thread_id, self.dead.shape[0], self.asp_solution.shape[0]
		for i in range(self.asp_solution.shape[0]):
			res = scipy.ndimage.interpolation.rotate(self.hrflat,-self.asp_solution[i,3],reshape=True,order=1,prefilter=False)
			coo = [self.asp_solution[i,1:3]]
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
			cent = np.zeros(2)
			cent[0] = np.floor(foc[0]-0.5)[1]
			cent[1] = np.floor(foc[0]-0.5)[0]
			#print cent
			half_len_x = int(res.shape[0]/2.)
			xl = cent[0]-half_len_x
			dxl = xl-0 if xl-0<=0 else 0
			xl = xl if xl>=0 else 0
			xh = cent[0]+(res.shape[0]-half_len_x)
			dxh = imsz[0]-xh if imsz[0]-xh<=0 else 0
			xh = xh if xh<imsz[0] else imsz[0]

			half_len_y = int(res.shape[1]/2.)
			yl = cent[1]-half_len_y
			dyl = yl-0 if yl-0<=0 else 0
			yl = yl if yl>=0 else 0
			yh = cent[1]+(res.shape[1]-half_len_y)
			dyh = imsz[1]-yh if imsz[1]-yh<=0 else 0
			yh = yh if yh<imsz[1] else imsz[1]
			#print xl,xh,yl,yh, res.shape
			#print dxl, dxh, dyl, dyh
			res_x, res_y = res.shape[0], res.shape[1]
			self.count[xl:xh, yl:yh] += res[-dxl:np.ceil(res_x+dxh), -dyl:np.ceil(res_y+dyh)]*0.005*(1-self.dead[i])
			if i%2000==0:
				print i
		print '%d done'%self.thread_id

class exp_dead_g_r(threading.Thread):
    def __init__(self, thread_id, imsz, wcs, hrflat, asp_solution, dead):
        threading.Thread.__init__(self)
        self.thread_id = thread_id
        self.count = np.zeros(imsz)
        self.wcs = wcs
        self.hrflat = hrflat
        self.asp_solution = asp_solution
        self.dead = dead
    def run(self):
		print '%d start'%self.thread_id, self.dead.shape[0], self.asp_solution.shape[0]

		hrflat = self.hrflat
		asp_solution = self.asp_solution
		count = self.count
		dead = self.dead
		rotate = scipy.ndimage.interpolation.rotate

		length = hrflat.shape[0]
		half_len = length/2.
		for i in range(asp_solution.shape[0]):
			res = rotate(hrflat,-asp_solution[i,3],reshape=True,order=1,prefilter=False)
			coo = [asp_solution[i,1:3]]
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
			cent = np.zeros(2)
			cent[0] = np.floor(foc[0]-0.5)[1]
			cent[1] = np.floor(foc[0]-0.5)[0]
			#print cent
			res_c_x = res.shape[0]/2.
			res_c_y = res.shape[1]/2.

			xl = cent[0]-half_len
			xh = cent[0]+half_len

			yl = cent[1]-half_len
			yh = cent[1]+half_len

			res_x, res_y = res.shape[0], res.shape[1]
			count[xl:xh, yl:yh] += res[res_c_x-half_len:res_c_x+half_len, res_c_y-half_len:res_c_y+half_len]*0.005*(1-dead[i])
			if i%2000==0:
				print i
		print '%d done'%self.thread_id

def exp_dead_new(file_num, imsz, wcs, hrflat, asp_solution, dead, return_dict):
	count = np.zeros(imsz)
	rotate = scipy.ndimage.interpolation.rotate

	length = hrflat.shape[0]
	half_len = length/2.
	print half_len
	coo = asp_solution[:,1:3]
	foc_list = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
	for i in range(asp_solution.shape[0]):
		res = rotate(hrflat,-asp_solution[i,3],reshape=False,order=1,prefilter=False)
		print res.shape
		foc = foc_list[i,:]#wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
		cent = np.zeros(2)
		cent[0] = np.floor(foc[1]-0.5)
		cent[1] = np.floor(foc[0]-0.5)
		#print cent
		res_c_x = res.shape[0]/2.
		res_c_y = res.shape[1]/2.

		xl = cent[0]-half_len
		xh = cent[0]+half_len

		yl = cent[1]-half_len
		yh = cent[1]+half_len

		res_x, res_y = res.shape[0], res.shape[1]
		#print cent[0], cent[1]
		#print xl,xh,yl,xh
		count[xl:xh, yl:yh] += res*0.005*(1-dead[i])#res[res_c_x-half_len:res_c_x+half_len, res_c_y-half_len:res_c_y+half_len]*0.005*(1-dead[i])
		if i%2000==0:
			print i
	print '%d done'%file_num
	return_dict[file_num] = count

if __name__ == '__main__':
	if False:
		scan_list = ['05','14','23','32','41','50','59','68']

		num_t = 7
		asp_solution_list = []
		ix_list = []
		dead_list = []
		for scan_num in scan_list:
			asp_solution_tmp = np.load('../data/photon_list/NUVPhoton%s_full_new_asp_cal_cata_wa.npy'%scan_num)
			asp_solution_list.append(asp_solution_tmp)			
			ix_tmp = (np.round(asp_solution_tmp[:,0]-asp_solution_tmp[0,0])+1).astype(int)
			dead_tmp = np.zeros(ix_tmp.shape[0])
			scst_tmp = pyfits.open('../data/scst/AIS_GAL_SCAN_000%s_0001-scst.fits'%scan_num)[1].data
			dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]
			dead_list.append(dead_tmp)


		asp_solution = np.concatenate(asp_solution_list, axis=0)
		dead = np.concatenate(dead_list, axis=0)

		asp_solutions = np.array_split(asp_solution, num_t)
		deads = np.array_split(dead, num_t)
		flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
		flat =flat_hdu[0].data
		size=flat_hdu[0].header['CDELT2']
		hrflat = scipy.ndimage.interpolation.zoom(flat,size/0.002,order=1,prefilter=False)

		hrflat = np.swapaxes(hrflat,0,1)

		skypos = [267.5, -26.0]
		skyrange = [25,20]
		tranges = []
		trange = [0, 0]
		imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
		wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)

		thread_list = []
		for i in range(num_t):
			thread = exp_dead(i, imsz, wcs, hrflat, asp_solutions[i], deads[i])
			thread.start()
			thread_list.append(thread)

		for t in thread_list:
			t.join()
		print 'all done'

		count = np.zeros(imsz)
		for t in thread_list:
			count += t.count

		tranges.append(trange)
		hdu = pyfits.PrimaryHDU(count)
		hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/count_map_05_68_gPr_flip_t_cut_new_dead.fits', clobber=False)

	if False:
		'''
		asp_solution1 = np.load('../data/photon_list/NUVPhoton05_full_new_asp_cal_cata.npy')
		asp_solution2 = np.load('../data/photon_list/NUVPhoton14_full_new_asp_cal_cata.npy')
		asp_solution = np.concatenate((asp_solution1,asp_solution2), axis=0)
		'''
		scan_list = ['05','14','23','32','41','50','59','68']
		num_t = 7
		asp_solution_list = []
		for scan_num in scan_list:
			asp_solution_list.append(np.load('../data/photon_list/NUVPhoton%s_full_new_asp_cal_cata_wa.npy'%scan_num))

		asp_solution = np.concatenate(asp_solution_list, axis=0)

		asp_solutions = np.array_split(asp_solution, num_t)

		asp_solution = asp_solutions[-1]

		flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
		flat =flat_hdu[0].data
		#flat = np.flipud(np.rot90(flat))
		size=flat_hdu[0].header['CDELT2']
		hrflat = scipy.ndimage.interpolation.zoom(flat,size/0.002,order=1,prefilter=False)

		hrflat = np.swapaxes(hrflat,0,1)

		skypos = [267.5, -26.0]
		skyrange = [25,20]
		tranges = []
		trange = [0, 0]
		imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
		count = np.zeros(imsz)
		wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)

		hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
		co_data = hdulist[1].data

		print asp_solution.shape

		for i in range(127952, asp_solution.shape[0]):
			res = scipy.ndimage.interpolation.rotate(hrflat,-asp_solution[i,3],reshape=True,order=1,prefilter=False)
			coo = [asp_solution[i,1:3]]
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
			cent = np.zeros(2)
			cent[0] = np.floor(foc[0]-0.5)[1]
			cent[1] = np.floor(foc[0]-0.5)[0]
			#print cent
			xl = cent[0]-res.shape[0]/2.
			dxl = xl-0 if xl-0<=0 else 0
			xl = xl if xl>=0 else 0
			xh = cent[0]+res.shape[0]/2.
			dxh = imsz[0]-xh if imsz[0]-xh<=0 else 0
			xh = xh if xh<imsz[0] else imsz[0]

			yl = cent[1]-res.shape[1]/2.
			dyl = yl-0 if yl-0<=0 else 0
			yl = yl if yl>=0 else 0
			yh = cent[1]+res.shape[1]/2.
			dyh = imsz[1]-yh if imsz[1]-yh<=0 else 0
			yh = yh if yh<imsz[1] else imsz[1]
			print xl,xh,yl,yh, res.shape
			res_x, res_y = res.shape[0], res.shape[1]
			print dxl, dxh, dyl, dyh
			count[xl:xh, yl:yh] += res[dxl:res_x+dxh, dyl:res_y+dyh]*0.005
			print i
			#if i>=140000:
			#	break

		tranges.append(trange)
		hdu = pyfits.PrimaryHDU(count)
		hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/count_map_05_14_gPr_res_r_flip.fits', clobber=False)

	if False:
		scan_list = ['05','14','23','32','41','50','59','68']
		num_t = 4
		asp_solution_list = []
		ix_list = []
		dead_list = []
		for scan_num in scan_list:
			asp_solution_tmp = np.load('../data/photon_list/NUVPhoton%s_full_new_asp_cal_cata_wa.npy'%scan_num)
			asp_solution_list.append(asp_solution_tmp)			
			ix_tmp = (np.round(asp_solution_tmp[:,0]-asp_solution_tmp[0,0])+1).astype(int)
			dead_tmp = np.zeros(ix_tmp.shape[0])
			scst_tmp = pyfits.open('../data/scst/AIS_GAL_SCAN_000%s_0001-scst.fits'%scan_num)[1].data
			dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]
			dead_list.append(dead_tmp)

		asp_solution = np.concatenate(asp_solution_list, axis=0)

		print asp_solution
		sky_data = SkyCoord(asp_solution[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
		print sky_data
		gal = sky_data.transform_to(Galactic)
		asp_solution[:,1:3] = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)

		print asp_solution
	
		dead = np.concatenate(dead_list, axis=0)

		asp_solutions = np.array_split(asp_solution, num_t)
		deads = np.array_split(dead, num_t)
		flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
		flat =flat_hdu[0].data
		size=flat_hdu[0].header['CDELT2']
		hrflat = scipy.ndimage.interpolation.zoom(flat,size/0.002,order=1,prefilter=False)

		hrflat = np.swapaxes(hrflat,0,1)

		hrflat = scipy.ndimage.interpolation.rotate(hrflat,59.,reshape=True,order=1,prefilter=False)

		skypos = [3.8, 0.]
		skyrange = [9., 24.]
		tranges = []
		trange = [0, 0]
		imsz = imagetools.deg2pix_g(skypos, skyrange, 0.002)
		wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)

		thread_list = []
		for i in range(num_t):
			thread = exp_dead_g_r(i, imsz, wcs, hrflat, asp_solutions[i], deads[i])
			thread.start()
			thread_list.append(thread)

		for t in thread_list:
			t.join()
		print 'all done'

		count = np.zeros(imsz)
		for t in thread_list:
			count += t.count

		tranges.append(trange)
		hdu = pyfits.PrimaryHDU(count)
		hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/count_map_05-68_gPr_flip_gal_dead_r_speed.fits', clobber=False)

#process
	if True:
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
		hrflat = scipy.ndimage.interpolation.zoom(flat,size/0.00166667,order=1,prefilter=False)
		print hrflat.shape
		hrflat = np.swapaxes(hrflat,0,1)
		print hrflat.shape
		hrflat = scipy.ndimage.interpolation.rotate(hrflat,59.,reshape=False,order=1,prefilter=False)
		print hrflat.shape
		skypos = [np.mean(gal_l), 0.]
		skyrange = [np.max(gal_l)-np.min(gal_l)+2., 24.]
		tranges = []
		trange = [0, 0]
		imsz = imagetools.deg2pix_g(skypos, skyrange, 0.00166667)
		wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.00166667)

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


