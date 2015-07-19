import numpy as np
from astropy.io import fits as pyfits
import imagetools
import pos_range
import scipy
import math
import threading

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
			res = scipy.ndimage.interpolation.rotate(hrflat,-self.asp_solution[i,3],reshape=True,order=1,prefilter=False)
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
			self.count[xl:xh, yl:yh] += res[dxl:res_x+dxh, dyl:res_y+dyh]*0.005
			if i%2000==0:
				print i
		print '%d done'%self.thread_id


if __name__ == '__main__':
	if True:
		scan_list = ['05','14','23','32','41','50','59','68']
		num_t = 7
		asp_solution_list = []
		for scan_num in scan_list:
			asp_solution_list.append(np.load('../data/photon_list/NUVPhoton%s_full_new_asp_cal_cata_wa.npy'%scan_num))

		asp_solution = np.concatenate(asp_solution_list, axis=0)

		asp_solutions = np.array_split(asp_solution, num_t)
		flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
		flat =flat_hdu[0].data
		size=flat_hdu[0].header['CDELT2']
		hrflat = scipy.ndimage.interpolation.zoom(flat,size/0.002,order=1,prefilter=False)

		hrflat = np.swapaxes(hrflat,0,1)

		skypos = [267.5, -29.0]
		skyrange = [26,24]
		tranges = []
		trange = [0, 0]
		imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
		wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)

		thread_list = []
		for i in range(num_t):
			thread = exp(i, imsz, wcs, hrflat, asp_solutions[i])
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
		hdulist.writeto('../fits/count_map_05_68_gPr_res_r_flip_t_cut_new.fits', clobber=False)

	if False:

		asp_solution1 = np.load('../data/photon_list/NUVPhoton05_full_new_asp_cal_cata.npy')
		asp_solution2 = np.load('../data/photon_list/NUVPhoton14_full_new_asp_cal_cata.npy')
		asp_solution = np.concatenate((asp_solution1,asp_solution2), axis=0)
		flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
		flat =flat_hdu[0].data
		#flat = np.flipud(np.rot90(flat))
		size=flat_hdu[0].header['CDELT2']
		hrflat = scipy.ndimage.interpolation.zoom(flat,size/0.002,order=1,prefilter=False)

		hrflat = np.swapaxes(hrflat,0,1)

		skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
		skypos = [268.5, -27.5]
		skyrange = [24,16]
		tranges = []
		trange = [0, 0]
		imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
		count = np.zeros(imsz)
		wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)

		hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
		co_data = hdulist[1].data

		print asp_solution.shape

		for i in range(asp_solution.shape[0]):
			res = scipy.ndimage.interpolation.rotate(hrflat,-asp_solution[i,3],reshape=True,order=1,prefilter=False)
			coo = [asp_solution[i,1:3]]
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
			cent = np.zeros(2)
			cent[0] = np.floor(foc[0]-0.5)[1]
			cent[1] = np.floor(foc[0]-0.5)[0]
			#print cent
			xl = cent[0]-res.shape[0]/2.
			xl = xl if xl>=0 else 0
			xh = cent[0]+res.shape[0]/2.
			xh = xh if xh<imsz[0] else imsz[0]

			yl = cent[1]-res.shape[1]/2.
			yl = yl if yl>=0 else 0
			yh = cent[1]+res.shape[1]/2.
			yh = yh if yh<imsz[1] else imsz[1]
			print xl,xh,yl,yh, res.shape
			count[xl:xh, yl:yh] += res*0.005
			print i
			#if i>=140000:
			#	break

		tranges.append(trange)
		hdu = pyfits.PrimaryHDU(count)
		hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
		hdulist = pyfits.HDUList([hdu])
		hdulist.writeto('../fits/count_map_05_14_gPr_res_r_flip.fits', clobber=False)


