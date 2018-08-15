import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import c3
from scipy import stats
import glob
import re
import os


def get_key(path):
	return int(re.split('_|/|\.', path)[-2])

def get_key_info(path):
	return int(re.split('_|\.', path)[-2])

def meta_plot(out, scan):
	#out='/scratch/dw1519/galex/co/co32_1-10_new'
	co_list = []
	info_list = []
	co_files =  sorted(glob.glob(out+"/[0-9]*.npy"), key=get_key)
	info_files = sorted(glob.glob(out+"/info_*"), key=get_key_info)
	print len(co_files)
	print len(info_files)
	'''
	for i in range(0, 27):
		co = np.load(out+'/radec/{0}.npy'.format((i+1)*100))/36000./800/0.001666*2400
		info = np.load(out+'/info_{0}.npy'.format((i+1)*100))[1:]
	'''
	for cf, inf in zip(co_files, info_files):
		co = np.load(cf)[1:]/36000./800/0.001666*2400
		info = np.load(inf)[1:]
		co_list.append(co)
		info_list.append(info)
	co = np.concatenate(co_list, axis=0)
	info = np.concatenate(info_list, axis=0)
	print co.shape
	print info.shape

	#ya hist
	mask = (co[:,1]>=-15) & (co[:,1]<=40)
	imsz = [32, 110]
	H,xedges,yedges=np.histogram2d(info[mask,1], co[mask,1],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 40, 0, 32], origin='lower', norm=LogNorm())

	'''
	plt.scatter(co[:,1], info[1:,1], s=0.5, alpha=0.005)
	plt.xlim(-15,40)
	plt.title('scan{0:>5} t={1}-{2}s'.format(239, i*50, (i+1)*50))
	plt.xlabel(r'$\Delta y$ [pixel]')
	plt.ylabel('YA')
	'''
	#ya = np.arange(30)+2
	#plt.plot(-0.0281724222131*ya+0.489247048825, ya, '-r') #239
	#plt.plot(-0.0426715749118*ya+0.714434337023, ya, '-r') #32
	#plt.title('scan{0:>5} t={1}-{2}s'.format(239, i*50, (i+1)*50))
	#plt.title(r'scan{0:>5} $\Delta y$=-0.043YA+0.714'.format(32))

	plt.xlim(-15,40)
	plt.ylim(0,32)
	plt.title('scan{0:>5}'.format(scan)) 
	plt.xlabel(r'$\Delta y$ [pixel]')
	plt.ylabel('YA')
	#plt.savefig(out+'/xya_{0}.png'.format((i+1)*100), dpi=190)
	plt.savefig(out+'/ya-y_tot.png', dpi=190)
	plt.clf()

	mask = (co[:,0]>=-15) & (co[:,0]<=15)
	imsz = [32, 60]
	H,xedges,yedges=np.histogram2d(info[mask,1], co[mask,0],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 15, 0, 32], origin='lower', norm=LogNorm())

	plt.xlim(-15,15)
	plt.ylim(0,32)
	plt.title('scan{0:>5}'.format(scan)) 
	plt.xlabel(r'$\Delta x$ [pixel]')
	plt.ylabel('YA')
	plt.savefig(out+'/ya-x_tot.png', dpi=190)
	plt.clf()

	#q hist
	mask = (co[:,0]>=-15) & (co[:,0]<=15) & (info[:,2]<24)
	imsz = [24, 60]
	H,xedges,yedges=np.histogram2d(info[mask,2], co[mask,0],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 15, 0, 24], origin='lower', norm=LogNorm())


	plt.xlim(-15,15)
	plt.ylim(0,24)
	plt.title('scan{0:>5}'.format(scan)) 
	plt.xlabel(r'$\Delta x$ [pixel]')
	plt.ylabel('q')
	plt.savefig(out+'/q-x_tot.png', dpi=190)
	plt.clf()

	mask = (co[:,1]>=-15) & (co[:,1]<=15) & (info[:,2]<24)
	imsz = [24, 60]
	H,xedges,yedges=np.histogram2d(info[mask,2], co[mask,1],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 15, 0, 24], origin='lower', norm=LogNorm())


	plt.xlim(-15,15)
	plt.ylim(0,24)
	plt.title('scan{0:>5}'.format(scan)) 
	plt.xlabel(r'$\Delta y$ [pixel]')
	plt.ylabel('q')
	plt.savefig(out+'/q-y_tot.png', dpi=190)
	plt.clf()

	#xa hist
	mask = (co[:,1]>=-15) & (co[:,1]<=15)
	imsz = [32, 60]
	H,xedges,yedges=np.histogram2d(info[mask,0], co[mask,1],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 15, 0, 32], origin='lower', norm=LogNorm())

	plt.xlim(-15,15)
	plt.ylim(0,32)
	plt.title('scan{0:>5}'.format(scan)) 
	plt.xlabel(r'$\Delta y$ [pixel]')
	plt.ylabel('XA')
	plt.savefig(out+'/xa-y_tot.png', dpi=190)
	plt.clf()

	mask = (co[:,0]>=-15) & (co[:,0]<=15)
	imsz = [32, 60]
	H,xedges,yedges=np.histogram2d(info[mask,0], co[mask,0],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 15, 0, 32], origin='lower', norm=LogNorm())

	plt.xlim(-15,15)
	plt.ylim(0,32)
	plt.title('scan{0:>5}'.format(scan)) 
	plt.xlabel(r'$\Delta x$ [pixel]')
	plt.ylabel('XA')
	plt.savefig(out+'/xa-x_tot.png', dpi=190)
	plt.clf()


def ya_fit(out, sacn):
	#out='/scratch/dw1519/galex/co/co32_1-10_new'
	co_list = []
	info_list = []
	co_files =  sorted(glob.glob(out+"/[0-9]*.npy"), key=get_key)
	info_files = sorted(glob.glob(out+"/info_*"), key=get_key_info)
	'''
	for i in range(0, 27):
		co = np.load(out+'/radec/{0}.npy'.format((i+1)*100))/36000./800/0.001666*2400
		info = np.load(out+'/info_{0}.npy'.format((i+1)*100))[1:]
	'''
	for cf, inf in zip(co_files, info_files):
		co = np.load(cf)[1:]/36000./800/0.001666*2400
		info = np.load(inf)[1:]
		co_list.append(co)
		info_list.append(info)

	co = np.concatenate(co_list, axis=0)
	info = np.concatenate(info_list, axis=0)
	mask1 = (co[:,1]>=-15) & (co[:,1]<=15)
	mask2 = (co[:,0]>=-15) & (co[:,0]<=15)
	centroids = []
	for ya in range(2,32):
		mask3 = info[:,1] == ya
		hist_size = [100,100]
		H,xedges,yedges=np.histogram2d(co[mask1&mask2&mask3,1], co[mask1&mask2&mask3,0],\
		                         bins=hist_size, range=([ [-15,15],[-15,15] ]))

		data = H.byteswap(True).newbyteorder()
		data = data.copy(order='C')
		data = data.byteswap(True).newbyteorder()
		t=50
		c_data = c3.find_centroid(data, t)
		if c_data is not None:
			cy, cx, max_value, flux = c_data
			centroids.append([-15.+cx*30./hist_size[0], -15.+cy*30./hist_size[0]])
	np.save(out+'/centroids.npy', centroids)
	ya = np.arange(30)+2
	centroids = np.array(centroids)
	plt.plot(ya, centroids[:,0], '.k')
	plt.xlabel('YA')
	plt.ylabel(r'mean $\Delta x$ [pixel]')
	plt.title(r'scan{0:>5}'.format(scan))
	plt.savefig(out+'/cx.png', dpi=190)
	plt.clf()
	plt.plot(ya, centroids[:,1], '.k')
	plt.xlabel('YA')
	plt.ylabel(r'mean $\Delta y$ [pixel]')
	plt.title(r'scan{0:>5}'.format(scan))
	plt.savefig(out+'/cy.png', dpi=190)
	plt.clf()

	c = np.load(out+'/centroids.npy')
	ya = np.arange(30)+2
	slope, intercept, r_value, p_value, std_err = stats.linregress(ya, c[:,1])
	print(slope, intercept, r_value, p_value, std_err)
	np.save(out+'/ya_fit.npy', [slope, intercept, r_value, p_value, std_err])

	#ya hist
	mask = (co[:,1]>=-15) & (co[:,1]<=40)
	imsz = [32, 110]
	H,xedges,yedges=np.histogram2d(info[mask,1], co[mask,1],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 40, 0, 32], origin='lower', norm=LogNorm())

	ya = np.arange(30)+2
	plt.plot(slope*ya+intercept, ya, '-r')
	plt.title(r'scan{0:>5} $\Delta y$={1:.3}YA+{2:.3}'.format(scan, slope, intercept))

	plt.xlim(-15,40)
	plt.ylim(0,32)
	plt.xlabel(r'$\Delta y$ [pixel]')
	plt.ylabel('YA')
	plt.savefig(out+'/ya-y_tot_fit.png', dpi=190)
	plt.clf()

def q_fit(out, sacn):
	#out='/scratch/dw1519/galex/co/co32_1-10_new'
	co_list = []
	info_list = []
	co_files =  sorted(glob.glob(out+"/[0-9]*.npy"), key=get_key)
	info_files = sorted(glob.glob(out+"/info_*"), key=get_key_info)
	'''
	for i in range(0, 27):
		co = np.load(out+'/radec/{0}.npy'.format((i+1)*100))/36000./800/0.001666*2400
		info = np.load(out+'/info_{0}.npy'.format((i+1)*100))[1:]
	'''
	for cf, inf in zip(co_files, info_files):
		co = np.load(cf)[1:]/36000./800/0.001666*2400
		info = np.load(inf)[1:]
		co_list.append(co)
		info_list.append(info)

	co = np.concatenate(co_list, axis=0)
	info = np.concatenate(info_list, axis=0)
	mask1 = (co[:,1]>=-15) & (co[:,1]<=15)
	mask2 = (co[:,0]>=-15) & (co[:,0]<=15)
	centroids = []
	for q in range(0,32):
		print('q={0}'.format(q))
		mask3 = info[:,2] == q
		hist_size = [100,100]
		H,xedges,yedges=np.histogram2d(co[mask1&mask2&mask3,1], co[mask1&mask2&mask3,0],\
		                         bins=hist_size, range=([ [-15,15],[-15,15] ]))

		data = H.byteswap(True).newbyteorder()
		data = data.copy(order='C')
		data = data.byteswap(True).newbyteorder()
		t=50
		c_data = c3.find_centroid(data, t)
		if c_data is not None:
			cy, cx, max_value, flux = c_data
			centroids.append([-15.+cx*30./hist_size[0], -15.+cy*30./hist_size[0]])
		else:
			centroids.append([0,0])
	np.save(out+'/centroids_q.npy', centroids)
	q = np.arange(32)
	centroids = np.array(centroids)
	plt.plot(q, centroids[:,0], '.k')
	plt.xlabel('q')
	plt.ylabel(r'mean $\Delta x$ [pixel]')
	plt.title(r'scan{0:>5}'.format(scan))
	plt.savefig(out+'/q_cx.png', dpi=190)
	plt.clf()
	plt.plot(q, centroids[:,1], '.k')
	plt.xlabel('q')
	plt.ylabel(r'mean $\Delta y$ [pixel]')
	plt.title(r'scan{0:>5}'.format(scan))
	plt.savefig(out+'/q_cy.png', dpi=190)
	plt.clf()

	'''
	c = np.load(out+'/centroids.npy')
	ya = np.arange(30)+2
	slope, intercept, r_value, p_value, std_err = stats.linregress(ya, c[:,1])
	print(slope, intercept, r_value, p_value, std_err)
	np.save(out+'/q_fit.npy', [slope, intercept, r_value, p_value, std_err])
	'''
	#q hist
	mask = (co[:,0]>=-15) & (co[:,0]<=15) & (info[:,2]<24)
	imsz = [24, 60]
	H,xedges,yedges=np.histogram2d(info[mask,2], co[mask,0],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 15, 0, 24], origin='lower', norm=LogNorm())

	q = np.arange(32)
	plt.plot(centroids[:,0], q, 'or')

	plt.xlim(-15,15)
	plt.ylim(0,24)
	plt.title('scan{0:>5}'.format(scan)) 
	plt.xlabel(r'$\Delta x$ [pixel]')
	plt.ylabel('q')
	plt.savefig(out+'/q-x_tot_fit.png', dpi=190)
	plt.clf()

	mask = (co[:,1]>=-15) & (co[:,1]<=15) & (info[:,2]<24)
	imsz = [24, 60]
	H,xedges,yedges=np.histogram2d(info[mask,2], co[mask,1],\
	                         bins=imsz)

	plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
				 extent=[-15, 15, 0, 24], origin='lower', norm=LogNorm())
	
	plt.plot(centroids[:,1], q, 'or')

	plt.xlim(-15,15)
	plt.ylim(0,24)
	plt.title('scan{0:>5}'.format(scan)) 
	plt.xlabel(r'$\Delta y$ [pixel]')
	plt.ylabel('q')
	plt.savefig(out+'/q-y_tot_fit.png', dpi=190)
	plt.clf()


if __name__ == '__main__':

	if False:
		void = ['co23_1-10', 'co32_1-10', 'co32-10', 'co32_1-10_new', 'co203_1-10', 'co239-10', 'co446_1-10', 
				'co464_1-10', 'co473_1-10', 'co806_1-10', 'co815_1-10', 'co1616_2-10', 'co1634_1-10', 'co1679_1-10', 
				'co2174_2-10', 'co2183_2-10', 'co2192_2-10', 'co2714_6-10', 'co2750_4-10', 'co3236_1-10', 
				'co3245_1-10', 'co3281_2-10']
		dirs = os.listdir('/scratch/dw1519/galex/co/')
		for d in dirs:
			if d not in void:
				print d
				out = '/scratch/dw1519/galex/co/'+d
				scan = re.split('co|_', d)[1]
				meta_plot(out, scan)
				ya_fit(out, scan)

	if False:
		void = ['co32_1-10', 'co32-10', 'co239-10']
		dirs = os.listdir('/scratch/dw1519/galex/co/')
		for d in dirs:
			if d not in void:
				print d
				out = '/scratch/dw1519/galex/co/'+d
				scan = re.split('co|_', d)[1]
				meta_plot(out, scan)
				ya_fit(out, scan)
				q_fit(out, scan)

	if False:
		dl = ['co203_1-10', 'co239_1-10', 'co1319_1-10', 'co1616_2-10', 'co1634_1-10']
		dl = ['co23_1-10', 'co32_1-10', 'co446_1-10','co464_1-10', 'co473_1-10', 'co815_1-10',\
			 'co1301_1-10', 'co2192_2-10']
		dl = ['co1679_1-10']
		dl = ['co806_1-10', 'co1310_1-10', 'co2174_2-10', 'co2183_2-10','co2714_6-10', 'co2750_4-10',\
			'co3236_1-10','co3245_1-10', 'co3281_2-10']
		for d in dl:
			out = '/scratch/dw1519/galex/co/'+d+'/test'
			scan = re.split('co|_', d)[1]
			meta_plot(out, scan)
			ya_fit(out, scan)
			q_fit(out, scan)



	if False:
		out='/home/dw1519/dw1519/galex/plots/co32-10/radec'
		co_list = []
		info_list = []
		for i in range(0, 28):
			co = np.load(out+'/{0}.npy'.format((i+1)*100))/36000./800/0.001666*2400
			info = np.load(out+'/info_{0}.npy'.format((i+1)*100))[1:]
			co_list.append(co)
			info_list.append(info)
		co = np.concatenate(co_list, axis=0)
		info = np.concatenate(info_list, axis=0)
		mask = (co[:,0]>=-15) & (co[:,0]<=15) & (info[:,2]<24)
		imsz = [24, 60]
		H,xedges,yedges=np.histogram2d(info[mask,2], co[mask,0],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-15, 15, 0, 24], origin='lower', norm=LogNorm())

		'''
		plt.scatter(co[:,1], info[1:,1], s=0.5, alpha=0.005)
		plt.xlim(-15,40)
		plt.title('scan{0:>5} t={1}-{2}s'.format(239, i*50, (i+1)*50))
		plt.xlabel(r'$\Delta y$ [pixel]')
		plt.ylabel('YA')
		'''
		#ya = np.arange(30)+2
		#plt.plot(-0.0281724222131*ya+0.489247048825, ya, '-r') #239
		#plt.plot(-0.0426715749118*ya+0.714434337023, ya, '-r') #32
		#plt.title('scan{0:>5} t={1}-{2}s'.format(239, i*50, (i+1)*50))
		#plt.title(r'scan{0:>5} $\Delta y$=-0.043YA+0.714'.format(32))

		plt.xlim(-15,15)
		plt.ylim(0,24)
		plt.title('scan{0:>5}'.format(32)) 
		plt.xlabel(r'$\Delta x$ [pixel]')
		plt.ylabel('XA')
		#plt.savefig(out+'/xya_{0}.png'.format((i+1)*100), dpi=190)
		plt.savefig(out+'/yq_tot.png', dpi=190)
		#print np.max(info[mask,2])


		plt.clf()
		plt.plot(H[0,:],'.c', label='q=0')
		plt.plot(H[1,:],'.y', label='q=1')
		plt.plot(H[2,:],'.g', label='q=2')
		plt.plot(H[3,:],'.b', label='q=3')
		plt.plot(H[4,:],'.k', label='q=4')
		plt.plot(H[5,:],'.r', label='q=5')
		plt.plot(H[6,:],'.m', label='q=6')
		plt.legend()
		plt.title('scan{0:>5}'.format(32)) 
		plt.xlabel(r'$\Delta x$ [pixel]')
		#plt.show()

		plt.savefig(out+'/yq_tot_hist.png', dpi=190)

	if False:
		out='/home/dw1519/dw1519/galex/plots/co239-10/radec'
		co_list = []
		info_list = []
		for i in range(0, 28):
			co = np.load(out+'/{0}.npy'.format((i+1)*100))/36000./800/0.001666*2400
			info = np.load(out+'/info_{0}.npy'.format((i+1)*100))[1:]
			co_list.append(co)
			info_list.append(info)
		co = np.concatenate(co_list, axis=0)
		info = np.concatenate(info_list, axis=0)
		plt.hist(info[:,1],32,alpha=0.5)
		plt.tile('{0}, photons'.format(32, np.sum(info[:,1]<2)/np.sum(info[:,1]>=0)))
		plt.savefig(out+'/ya_tot_hist.png', dpi=190)

	if False:
		scan = '2813'
		name = 'AIS_GAL_SCAN_0'+scan+'_0003'
		out = '../data/'+name+'-cal-sec'
		data = np.load('../data/'+name+'-cal-sec/photon_match_d_081017.npy')
		print scan
		data_mask = data[:,-3]>-1
		data = data[data_mask]
		dy = data[:,-1]/36000./800/0.001666*2400
		dx = data[:,-2]/36000./800/0.001666*2400
		ya = data[:,-5]
		q = data[:,-6]
		ya = data[:,-7]
		xa = data[:,-8]
		print dy
		print np.max(dy)
		print np.min(dy)
		print np.percentile(dy, 80)
		mask = (dy>=-8) & (dy<=8) & (q<24)
		imsz = [24, 64]
		H,xedges,yedges=np.histogram2d(q[mask], dy[mask],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-8, 8, 0, 24], origin='lower', norm=LogNorm())

		dy_cs = []
		for i in range(24):
			q_mask = q[mask] == i
			dy_c = np.median(dy[mask][q_mask])
			dy_cs.append(dy_c)
		plt.plot(dy_cs, np.arange(24), '.r')

		plt.xlim(-8,8)
		plt.ylim(0,24)
		plt.title('scan{0:>5}'.format(scan)) 
		plt.xlabel(r'$\Delta y$ [pixel]')
		plt.ylabel('q')
		plt.savefig(out+'/q-y_tot.png', dpi=190)
		plt.clf()
		dy_q = np.array([np.arange(24),dy_cs]).T
		np.save(out+'/q-y.npy', dy_q)

		plt.plot(np.arange(24), dy_cs, '.r')
		plt.savefig(out+'/q-y_c.png', dpi=190)
		plt.clf()

		mask = (dx>=-8) & (dx<=8) & (q<24)
		imsz = [24, 64]
		H,xedges,yedges=np.histogram2d(q[mask], dx[mask],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-8, 8, 0, 24], origin='lower', norm=LogNorm())


		dy_cs = []
		for i in range(24):
			q_mask = q[mask] == i
			dy_c = np.median(dx[mask][q_mask])
			dy_cs.append(dy_c)
		plt.plot(dy_cs, np.arange(24), '.r')

		plt.xlim(-8,8)
		plt.ylim(0,24)
		plt.title('scan{0:>5}'.format(scan)) 
		plt.xlabel(r'$\Delta x$ [pixel]')
		plt.ylabel('q')
		plt.savefig(out+'/q-x_tot.png', dpi=190)
		plt.clf()
		dx_q = np.array([np.arange(24),dy_cs]).T
		np.save(out+'/q-x.npy', dx_q)

		plt.plot(np.arange(24), dy_cs, '.r')
		plt.savefig(out+'/q-x_c.png', dpi=190)
		plt.clf()

		q_mask = q<24
		dy[q_mask] = dy[q_mask]-dy_q[data[q_mask,-6].astype(int),-1]
		dx[q_mask] = dx[q_mask]-dx_q[data[q_mask,-6].astype(int),-1]

		mask = (dy>=-8) & (dy<=8)
		imsz = [30, 64]
		H,xedges,yedges=np.histogram2d(ya[mask], dy[mask],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-8, 8, 2, 32], origin='lower', norm=LogNorm())

		dy_cs = []
		for i in range(2,32):
			ya_mask = ya[mask] == i
			dy_c = np.median(dy[mask][ya_mask])
			dy_cs.append(dy_c)
		plt.plot(dy_cs, np.arange(30)+2, '.r')

		plt.xlim(-8,8)
		plt.ylim(2,32)
		plt.title('scan{0:>5}'.format(scan)) 
		plt.xlabel(r'$\Delta y$ [pixel]')
		plt.ylabel('ya')
		plt.savefig(out+'/ya-y_tot.png', dpi=190)
		plt.clf()

		plt.plot(np.arange(30)+2, dy_cs, '.k')
		plt.savefig(out+'/ya-y_c.png', dpi=190)
		plt.clf()
		np.save(out+'/ya-y.npy', np.array([np.arange(30)+2,dy_cs]).T)

		mask = (dy>=-8) & (dy<=8)
		imsz = [30, 64]
		H,xedges,yedges=np.histogram2d(ya[mask], dx[mask],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-8, 8, 2, 32], origin='lower', norm=LogNorm())

		dy_cs = []
		for i in range(2,32):
			ya_mask = ya[mask] == i
			dy_c = np.median(dx[mask][ya_mask])
			dy_cs.append(dy_c)
		plt.plot(dy_cs, np.arange(30)+2, '.r')

		plt.xlim(-8,8)
		plt.ylim(2,32)
		plt.title('scan{0:>5}'.format(scan)) 
		plt.xlabel(r'$\Delta x$ [pixel]')
		plt.ylabel('ya')
		plt.savefig(out+'/ya-x_tot.png', dpi=190)
		plt.clf()

		plt.plot(np.arange(30)+2, dy_cs, '.k')
		plt.savefig(out+'/ya-x_c.png', dpi=190)
		plt.clf()
		np.save(out+'/ya-x.npy', np.array([np.arange(30)+2,dy_cs]).T)


	if False:
		'''
		scan = '2561'
		name = 'AIS_GAL_SCAN_0'+scan+'_0002'
		out = '../data/'+name+'-cal-sec'
		data = np.load('../data/'+name+'-cal-sec/photon_match_d_081017.npy')
		print scan
		data_mask = data[:,-3]>-1
		data = data[data_mask]
		'''
		dy = data[:,-1]/36000./800/0.001666*2400
		dx = data[:,-2]/36000./800/0.001666*2400
		q = data[:,-6].astype(int)
		ya = data[:,-7].astype(int)
		xa = data[:,-8].astype(int)

		dy_q = np.load(out+'/q-y.npy')
		dx_q = np.load(out+'/q-x.npy')
		q_mask = q<24
		dy[q_mask] = dy[q_mask]-dy_q[data[q_mask,-6].astype(int),-1]
		dx[q_mask] = dx[q_mask]-dx_q[data[q_mask,-6].astype(int),-1]

		dy_ya = np.load(out+'/ya-y.npy')
		dx_ya = np.load(out+'/ya-x.npy')
		dy = dy-dy_ya[data[:,-7].astype(int)-2, -1]
		dx = dx-dx_ya[data[:,-7].astype(int)-2, -1]
		print dy
		print np.max(dy)
		print np.min(dy)
		print np.percentile(dy, 80)
		mask = (dy>=-8) & (dy<=8) & (q<24)
		imsz = [24, 64]
		H,xedges,yedges=np.histogram2d(q[mask], dy[mask],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-8, 8, 0, 24], origin='lower', norm=LogNorm())

		dy_cs = []
		for i in range(24):
			q_mask = q[mask] == i
			dy_c = np.median(dy[mask][q_mask])
			dy_cs.append(dy_c)
		plt.plot(dy_cs, np.arange(24), '.r')

		plt.xlim(-8,8)
		plt.ylim(0,24)
		plt.title('scan{0:>5}'.format(scan)) 
		plt.xlabel(r'$\Delta y$ [pixel]')
		plt.ylabel('q')
		plt.savefig(out+'/q-y_tot_new.png', dpi=190)
		plt.clf()
		#np.save(out+'/q-y.npy', np.array([np.arange(24),dy_cs]).T)

		plt.plot(np.arange(24), dy_cs, '.r')
		plt.savefig(out+'/q-y_c_new.png', dpi=190)
		plt.clf()

		mask = (dx>=-8) & (dx<=8) & (q<24)
		imsz = [24, 64]
		H,xedges,yedges=np.histogram2d(q[mask], dx[mask],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-8, 8, 0, 24], origin='lower', norm=LogNorm())


		dy_cs = []
		for i in range(24):
			q_mask = q[mask] == i
			dy_c = np.median(dx[mask][q_mask])
			dy_cs.append(dy_c)
		plt.plot(dy_cs, np.arange(24), '.r')

		plt.xlim(-8,8)
		plt.ylim(0,24)
		plt.title('scan{0:>5}'.format(scan)) 
		plt.xlabel(r'$\Delta x$ [pixel]')
		plt.ylabel('q')
		plt.savefig(out+'/q-x_tot_new.png', dpi=190)
		plt.clf()
		#np.save(out+'/q-x.npy', np.array([np.arange(24),dy_cs]).T)

		plt.plot(np.arange(24), dy_cs, '.r')
		plt.savefig(out+'/q-x_c_new.png', dpi=190)
		plt.clf()


		mask = (dy>=-8) & (dy<=8)
		imsz = [30, 64]
		H,xedges,yedges=np.histogram2d(ya[mask], dy[mask],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-8, 8, 2, 32], origin='lower', norm=LogNorm())

		dy_cs = []
		for i in range(2,32):
			ya_mask = ya[mask] == i
			dy_c = np.median(dy[mask][ya_mask])
			dy_cs.append(dy_c)
		plt.plot(dy_cs, np.arange(30)+2, '.r')

		plt.xlim(-8,8)
		plt.ylim(2,32)
		plt.title('scan{0:>5}'.format(scan)) 
		plt.xlabel(r'$\Delta y$ [pixel]')
		plt.ylabel('ya')
		plt.savefig(out+'/ya-y_tot_new.png', dpi=190)
		plt.clf()

		plt.plot(np.arange(30)+2, dy_cs, '.k')
		plt.savefig(out+'/ya-y_c_new.png', dpi=190)
		plt.clf()
		#np.save(out+'/ya-y.npy', np.array([np.arange(30)+2,dy_cs]).T)

		mask = (dy>=-8) & (dy<=8)
		imsz = [30, 64]
		H,xedges,yedges=np.histogram2d(ya[mask], dx[mask],\
		                         bins=imsz)

		plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
					 extent=[-8, 8, 2, 32], origin='lower', norm=LogNorm())

		dy_cs = []
		for i in range(2,32):
			ya_mask = ya[mask] == i
			dy_c = np.median(dx[mask][ya_mask])
			dy_cs.append(dy_c)
		plt.plot(dy_cs, np.arange(30)+2, '.r')

		plt.xlim(-8,8)
		plt.ylim(2,32)
		plt.title('scan{0:>5}'.format(scan)) 
		plt.xlabel(r'$\Delta x$ [pixel]')
		plt.ylabel('ya')
		plt.savefig(out+'/ya-x_tot_new.png', dpi=190)
		plt.clf()

		plt.plot(np.arange(30)+2, dy_cs, '.k')
		plt.savefig(out+'/ya-x_c_new.png', dpi=190)
		plt.clf()
		#np.save(out+'/ya-x.npy', np.array([np.arange(30)+2,dy_cs]).T)

	if True:
		scans = ['2552', '2561', '2570', '2615', '2642', '2921', '3497']
		scans = ['2813']
		for scan in scans:
			#scan = '2552'
			name = 'AIS_GAL_SCAN_0'+scan+'_0003'
			out = '../data/'+name+'-cal-sec'

			try:
				dy_q = np.load(out+'/q-y.npy')
				dx_q = np.load(out+'/q-x.npy')

				dy_ya = np.load(out+'/ya-y.npy')
				dx_ya = np.load(out+'/ya-x.npy')
			except IOError:
				continue

			np.save(out+'/q_corr.npy', np.column_stack([dx_q,dy_q[:,-1]]))
			np.save(out+'/ya_corr.npy', np.column_stack([dx_ya,dy_ya[:,-1]]))


	if False:
		scans = ['2552', '2561', '2570', '2615', '2642', '2921', '3497']
		fig, ax = plt.subplots(2,2)
		for scan in scans:
			#scan = '2552'
			for num in range(1,4):
				name = 'AIS_GAL_SCAN_0'+scan+'_000{0}'.format(num)
				out = '../data/'+name+'-cal-sec'

				try:
					q_corr = np.load(out+'/q_corr.npy')
					ya_corr = np.load(out+'/ya_corr.npy')
				except IOError:
					continue

				ax[0,0].plot(q_corr[:,0], q_corr[:,1], label=scan+'_000{0}'.format(num))
				ax[0,0].set_xlabel('q')
				ax[0,0].set_ylabel('dx')

				ax[0,1].plot(q_corr[:,0], q_corr[:,2], label=scan+'_000{0}'.format(num))
				ax[0,1].set_xlabel('q')
				ax[0,1].set_ylabel('dy')

				ax[1,0].plot(ya_corr[:,0], ya_corr[:,1], label=scan+'_000{0}'.format(num))
				ax[1,0].set_xlabel('ya')
				ax[1,0].set_ylabel('dx')

				ax[1,1].plot(ya_corr[:,0], ya_corr[:,2], label=scan+'_000{0}'.format(num))
				ax[1,1].set_xlabel('ya')
				ax[1,1].set_ylabel('dy')

		ax[0,0].legend()
		ax[0,1].legend()
		ax[1,0].legend()
		ax[1,1].legend()

		plt.savefig('../data/compare.png', dpi=190)


	if False:
		scans = ['2552', '2561', '2570', '2615', '2642', '2921', '3497']
		fig, ax = plt.subplots(2,1)
		fig.set_size_inches(5, 8)
		for scan in scans:
			#scan = '2552'
			for num in range(1,4):
				name = 'AIS_GAL_SCAN_0'+scan+'_000{0}'.format(num)
				out = '../data/'+name+'-cal-sec'

				try:
					q_corr = np.load(out+'/q_corr.npy')
					#ya_corr = np.load(out+'/ya_corr.npy')
				except IOError:
					continue

				ax[0].plot(q_corr[:,0], q_corr[:,1], label=scan+'_000{0}'.format(num))
				ax[0].set_xlabel('q')
				ax[0].set_ylabel('dx')

				ax[1].plot(q_corr[:,0], q_corr[:,2], label=scan+'_000{0}'.format(num))
				ax[1].set_xlabel('q')
				ax[1].set_ylabel('dy')


		#ax[0,0].legend()
		#ax[0,1].legend()

		plt.savefig('../data/q_compare.png', dpi=190)


