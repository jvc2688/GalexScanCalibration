from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import load_data
from matplotlib.colors import LogNorm
#plt.style.use('seaborn-paper')

def plot(date, scan):
	#date = '08-21-2017'
	print(scan)
	tycho2 = pyfits.open('/scratch/dw1519/galex/data/tycho2.fits')[1].data
	cfile='/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
	df = load_data.load_catalog(cfile)
	c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
	'''
	if scan=='3488':
		cut = df['gb']>1.
		c = c[cut]
		df = df[cut]
	'''
	print(c.shape)
	catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
	idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
	mask=d2d<0.001*u.degree
	N = np.sum(mask)
	print(N)

	f = plt.figure(figsize=(12,10))
	ax = plt.subplot(2, 2, 1)
	ax.grid(False)
	plt.scatter((c[mask].l.deg-catalog[idx[mask]].l.deg)*3600,\
		(c[mask].b.deg-catalog[idx[mask]].b.deg)*3600,\
		 c=df['nuv'][mask], alpha=0.2, edgecolors='face', cmap=plt.get_cmap('jet'))
	plt.axis('equal')
	plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
	plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
	#plt.title('0023 NUV<15, star={0} {1:.0f}%, photon={2} {3:.0f}%'.format(len(bright), len(bright)/total_star*100,\
	#												 len(bright_co), len(bright_co)/total_photon*100))
	plt.xlim(-3.0,3.0)
	plt.ylim(-3.0,3.0)
	plt.ylabel('gb(SEx-T2)')
	plt.title('gl(SEx-T2), N={0}'.format(N))
	plt.colorbar()


	ais = pyfits.open('/scratch/dw1519/galex/data/MyTable_4_jvc2688.fit')[1].data
	ais_cat = SkyCoord(ais['glon']*u.degree, ais['glat']*u.degree, frame='galactic')
	idx, d2d, d3d = c.match_to_catalog_sky(ais_cat)  
	mask=d2d<0.001*u.degree
	N = np.sum(mask)
	print('ais match:{0}'.format(N))
	plt.subplot(2, 2, 2)
	plt.scatter(ais['nuv'][idx[mask]],\
		df['nuv'][mask]-ais['nuv'][idx[mask]],\
		 c=ais['glat'][idx[mask]], alpha=0.5, edgecolors='face', s=0.2, cmap=plt.get_cmap('jet'), vmin=-10, vmax=10)
	plt.axhline(y=0, color='k', linestyle='dashed', linewidth=1)
	nuv_list = [[14,15],[15,16],[16,17],[17,18],[18,19],[19,20]]
	mask_list = [(e[0]<=ais['nuv'][idx[mask]]) & (e[1]>ais['nuv'][idx[mask]]) for e in nuv_list]
	x = [np.median(ais['nuv'][idx[mask]][m]) for m in mask_list]
	y = [np.median((df['nuv'][mask]-ais['nuv'][idx[mask]])[m]) for m in mask_list]
	ye = [np.std((df['nuv'][mask]-ais['nuv'][idx[mask]])[m]) for m in mask_list]
	plt.errorbar(x, y, yerr=ye)
	plt.xlim(12,22)
	plt.ylabel('NUV (SEx-AIS)')
	plt.title('NUV AIS, N={0}'.format(N))
	plt.ylim(-2,2)
	plt.colorbar()

	plt.subplot(2, 2, 3)
	plt.hist(df['nuv'],bins=18,range=(12,22))
	plt.title('NUV SExtractor')
	plt.xlim(12,22)

	plt.subplot(2, 2, 4)
	plt.plot(df['nuv'], df['FWHM_IMAGE'],'.k',markersize=0.7)
	plt.title('NUV SExtractor, N={0}'.format(len(df)))
	plt.ylabel('FWHM [pix]')
	plt.ylim(0,8)
	plt.xlim(12,22)

	plt.suptitle(scan)
	plt.tight_layout()
	#plt.subplots_adjust(top=0.9)

	plt.savefig('/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/'+scan+'_match_new.png', dpi=190)
	plt.clf()



if __name__ == '__main__':
	if True:
		scans = ['0023', '0032', '0203', '0239', '0446', '0464', '0473', '0806', '0815', '1301', '1310', '1319',\
				'1616', '1679', '2192', '2714', '2750', '3236', '3245', '3281']

		scans = ['0050']#['0050']#['0-10ds'] #['0050']
		#scans = ['1634', '2183', '2552', '2561', '2570', '2615', '2642', '2921', '3497']
		#scans = ['2174']
		date = '11-07-2017'
		for scan in scans:
			print(scan)
			tycho2 = pyfits.open('../data/tycho2.fits')[1].data
			cfile='/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm_in.txt'
			df = load_data.load_catalog(cfile)
			c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
			print(c.shape)
			catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
			idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
			mask=d2d<0.001*u.degree
			N = np.sum(mask)
			print(N)

			f = plt.figure(figsize=(8,5))



			ais = pyfits.open('../data/MyTable_4_jvc2688.fit')[1].data
			ais_cat = SkyCoord(ais['glon']*u.degree, ais['glat']*u.degree, frame='galactic')
			idx, d2d, d3d = c.match_to_catalog_sky(ais_cat)  
			mask=(d2d<0.001*u.degree) #& (df['gl']<50)
			N = np.sum(mask)
			print('ais match:{0}'.format(N))
			#plt.subplot(2, 2, 2)
			gl = df['gl'][mask]
			gb = df['gb'][mask]
			plt.scatter(gl,\
				df['nuv'][mask]-ais['nuv'][idx[mask]],\
				alpha=0.8, edgecolors='face', s=0.5, cmap=plt.get_cmap('jet'), vmin=-10, vmax=10)
				# c=ais['glat'][idx[mask]], alpha=0.8, edgecolors='face', s=0.5, cmap=plt.get_cmap('jet'), vmin=-10, vmax=10)
			plt.axhline(y=0, color='k', linestyle='dashed', linewidth=1)
			#nuv_list = [[14,15],[15,16],[16,17],[17,18],[18,19],[19,20]]
			#mask_list = [(e[0]<=ais['nuv'][idx[mask]]) & (e[1]>ais['nuv'][idx[mask]]) for e in nuv_list]
			#x = [np.median(ais['nuv'][idx[mask]][m]) for m in mask_list]
			#y = [np.median((df['nuv'][mask]-ais['nuv'][idx[mask]])[m]) for m in mask_list]
			#ye = [np.std((df['nuv'][mask]-ais['nuv'][idx[mask]])[m]) for m in mask_list]
			#plt.errorbar(x, y, yerr=ye)
			#plt.xlim(0,4)
			plt.ylabel('NUV (SEx-AIS)')
			plt.xlabel('gl')
			plt.title('NUV AIS, N={0}'.format(N))
			plt.ylim(-2,2)
			#plt.colorbar()


			#plt.suptitle(scan)
			plt.tight_layout()
			#plt.subplots_adjust(top=0.9)

			plt.savefig('/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/'+scan+'_match_exp.png', dpi=190)
			plt.clf()
			np.save('/scratch/dw1519/galex/data/star_photon/iterate3/stars_in.npy', 
				np.column_stack([gl,df['nuv'][mask]-ais['nuv'][idx[mask]], ais['nuv'][idx[mask]], df['FWHM_IMAGE'][mask]]))
			np.save('/scratch/dw1519/galex/data/star_photon/iterate3/stars_in_all.npy',
				np.column_stack([df['gl'], df['nuv'], df['FWHM_IMAGE']]))





