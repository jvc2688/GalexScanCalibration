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
	if False:
		name_list = ['0023', '0032', '0203', '0239', '0446', '0464', '0473', '0806', '0815', '1301', '1310', '1319',\
					'1616', '1634', '1679', '2174', '2183', '2192', '2714', '2750', '3236', '3245', '3281']
		name_list = ['2552', '2561', '2570', '2615', '2642', '2921', '3497']
		tycho2 = pyfits.open('../data/tycho2.fits')[1].data
		for scan in name_list:
			#scan = '0023'
			print(scan)
			cfile='/scratch/dw1519/galex/fits/scan_map/catalog/08-10-2017_1/starcat_{0}mapweight_fec_fwhm.txt'.format(scan)
			try:
				df = load_data.load_catalog(cfile)
			except IOError:
				print('skip')
				continue
			c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
			catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
			idx, d2d, d3d = c.match_to_catalog_sky(catalog)
			mask=d2d<0.001*u.degree
			print(np.sum(mask))
			dtype = np.dtype([('tycho_num', int), ('Glon', '>f4'), ('Glat', '>f4'), ('RAJ2000', '>f4'), ('DEJ2000', '>f4'),
							 ('flux', float), ('nuv', float), ('gl', float), ('gb', float)])
			matched_tycho = tycho2[idx[mask]]
			matched_df = df[mask]
			matched_catalog = np.core.records.fromarrays(np.array([idx[mask], matched_tycho['Glon'], matched_tycho['Glat'],
							matched_tycho['RAJ2000'], matched_tycho['DEJ2000'], np.array(matched_df['FLUX_AUTO']),
							np.array(matched_df['nuv']), np.array(matched_df['gl']),
							np.array(matched_df['gb'])]), dtype=dtype)
			print(matched_catalog.shape)
			np.save('/scratch/dw1519/galex/fits/scan_map/catalog/08-10-2017_1/starcat_{0}mapweight_fec_fwhm.npy'.format(scan),
					matched_catalog)


	if False:
		tycho2 = pyfits.open('../data/tycho2.fits')[1].data
		cfile='/scratch/dw1519/galex/data/extraction/Starcat-07-22-2017/starcat_0023mapweight_fec_fwhm.txt'
		df = load_data.load_catalog(cfile)
		c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
		catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
		idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
		mask=d2d<0.001*u.degree

		co = np.load('../data/AIS_GAL_SCAN_00023_0001-cal-sec/co_match_072217.npy')
		photon_nuv = np.zeros(co.shape[0])
		nuv = np.array(df['nuv'][mask])
		for i in range(idx[mask].shape[0]):
			photon_nuv[co[:,-1]==idx[mask][i]] = nuv[i]

		if True:
			total_star = float(np.sum(mask))
			total_photon = float(co.shape[0])

			bright_mask = df['nuv'][mask]<14.4
			bright = set(idx[mask][bright_mask])
			bright_co = co[photon_nuv<14.4]#np.array([a for a in co if a[-1] in bright])
			weights = np.power(10., -photon_nuv[photon_nuv<14.4]/2.5)
			load_data.plot_co(bright_co, 0.00008, [0,0], [0.0045,0.0045], weights)
			plt.scatter((c[mask][bright_mask].l.deg-catalog[idx[mask]][bright_mask].l.deg)*3600,\
				(c[mask][bright_mask].b.deg-catalog[idx[mask]][bright_mask].b.deg)*3600,\
				 c=df['nuv'][mask][bright_mask], alpha=0.2, edgecolors='face')

			plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
			plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
			plt.title('0023 NUV<15, star={0} {1:.0f}%, photon={2} {3:.0f}%'.format(len(bright), len(bright)/total_star*100,\
															 len(bright_co), len(bright_co)/total_photon*100))
			plt.xlim(-4.5,4.5)
			plt.ylim(-4.5,4.5)
			plt.colorbar()
			plt.savefig('../data/AIS_GAL_SCAN_00023_0001-cal-sec/bright.png', dpi=190)
			plt.clf()

			'''
			median_mask = (df['nuv'][mask]>=14.4) & (df['nuv'][mask]<15.6)
			median = set(idx[mask][median_mask])
			median_co = co[(photon_nuv>=14.4) & (photon_nuv<15.6)]#np.array([a for a in co if a[-1] in median])
			load_data.plot_co(median_co, 0.00008, [0,0], [0.0045,0.0045])

			plt.scatter((c[mask][median_mask].l.deg-catalog[idx[mask]][median_mask].l.deg)*3600,\
				(c[mask][median_mask].b.deg-catalog[idx[mask]][median_mask].b.deg)*3600,\
				 c=df['nuv'][mask][median_mask], alpha=0.2, edgecolors='face')
			plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
			plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
			plt.title('0023 14<=NUV<=15.6, star={0} {1:.0f}%, photon={2} {3:.0f}%'.format(len(median), len(median)/total_star*100,\
															 len(median_co), len(median_co)/total_photon*100))		
			plt.xlim(-4.5,4.5)
			plt.ylim(-4.5,4.5)
			plt.colorbar()
			plt.show()
			'''

			faint_mask = (df['nuv'][mask]>=14.4) #& (df['nuv'][mask]<15.6)
			faint = set(idx[mask][faint_mask])
			faint_co = co[photon_nuv>=14.4]#np.array([a for a in co if a[-1] in faint])
			weights = np.power(10., photon_nuv[photon_nuv>=14.4]/2.5)/100000.

			load_data.plot_co(faint_co, 0.00008, [0,0], [0.0045,0.0045], weights)

			plt.scatter((c[mask][faint_mask].l.deg-catalog[idx[mask]][faint_mask].l.deg)*3600,\
				(c[mask][faint_mask].b.deg-catalog[idx[mask]][faint_mask].b.deg)*3600,\
				 c=df['nuv'][mask][faint_mask], alpha=0.2, edgecolors='face')
			plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
			plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
			plt.title('0023 NUV>=15, star={0} {1:.0f}%, photon={2} {3:.0f}%'.format(len(faint), len(faint)/total_star*100,\
															 len(faint_co), len(faint_co)/total_photon*100))
			plt.xlim(-4.5,4.5)
			plt.ylim(-4.5,4.5)
			plt.colorbar()
			plt.savefig('../data/AIS_GAL_SCAN_00023_0001-cal-sec/faint.png', dpi=190)
			plt.clf()

		if False:
			imsz = [26, 20]
			H,xedges,yedges=np.histogram2d(co[:,4], photon_nuv, bins=imsz)
			print H.shape
			H = H/(np.sum(H,axis=0).astype(float))
			plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='auto',\
						 extent=[10, 20, 6, 32], origin='lower')

			plt.xlabel('NUV')
			plt.ylabel('q')
			plt.title('0023')
			plt.savefig('../data/AIS_GAL_SCAN_00023_0001-cal-sec/q-nuv.png', dpi=190)
			plt.clf()

	if True:
		scans = ['0023', '0032', '0203', '0239', '0446', '0464', '0473', '0806', '0815', '1301', '1310', '1319',\
				'1616', '1679', '2192', '2714', '2750', '3236', '3245', '3281']

		scans = ['0050']#['2561']
		#scans = ['1634', '2183', '2552', '2561', '2570', '2615', '2642', '2921', '3497']
		#scans = ['2174']
		date = '08-21-2017'#'08-21-2017'
		for scan in scans:
			print(scan)
			tycho2 = pyfits.open('../data/tycho2.fits')[1].data
			cfile='/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
			df = load_data.load_catalog(cfile)
			c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
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


			ais = pyfits.open('../data/MyTable_4_jvc2688.fit')[1].data
			ais_cat = SkyCoord(ais['glon']*u.degree, ais['glat']*u.degree, frame='galactic')
			idx, d2d, d3d = c.match_to_catalog_sky(ais_cat)  
			mask=d2d<0.001*u.degree
			N = np.sum(mask)
			print('ais match:{0}'.format(N))
			plt.subplot(2, 2, 2)
			plt.scatter(ais['nuv'][idx[mask]],\
				df['nuv'][mask]-ais['nuv'][idx[mask]],\
				 c=ais['glon'][idx[mask]], alpha=0.5, edgecolors='face', s=0.2, cmap=plt.get_cmap('jet'))
				 #c=ais['glat'][idx[mask]], alpha=0.5, edgecolors='face', s=0.2, cmap=plt.get_cmap('jet'), vmin=-10, vmax=10)
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

			plt.savefig('/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/'+scan+'_match_new_gl.png', dpi=190)
			plt.clf()

			'''
			cfile='/scratch/dw1519/galex/data/extraction/Starcat-07-22-2017/starcat_'+scan+'mapweight_fec_fwhm.txt'
			df = load_data.load_catalog(cfile)
			c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
			print(c.shape)
			catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
			idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
			mask=d2d<0.001*u.degree
			N = np.sum(mask)
			print(np.sum(mask))
			f = plt.figure(figsize=(12,10))
			plt.subplot(2, 2, 1)
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

			#ais = pyfits.open('../data/MyTable_2_jvc2688.fit')[1].data
			#ais_cat = SkyCoord(ais['glon']*u.degree, ais['glat']*u.degree, frame='galactic')
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
			plt.savefig('/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/'+scan+'_match_old.png', dpi=190)
			plt.clf()
			'''




