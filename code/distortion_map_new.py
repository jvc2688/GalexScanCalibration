import matplotlib.pyplot as plt
import numpy as np
import c3

def split_seq(seq, size):
	newseq = []
	splitsize = 1.0/size*len(seq)
	for i in range(size):
		newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
	return newseq


if __name__ == '__main__':
	suffix = 'xa_2'
	scan_list = []
	with open('../check/distortion_list', 'r') as f:
		scan_list = f.read().splitlines()
	scan_list = sorted(scan_list)[0:30]
	print len(scan_list)

	split = split_seq(scan_list, 1)
	print len(split)

	hist_size = [100,100]
	detector_size = [50,50]

	for name_list in split:
		print name_list
		hist_list = []
		for i in range(detector_size[0]):
			for j in range(detector_size[1]):
				hist_list.append(np.zeros(hist_size))

		for name in name_list: 
			hists = np.load('../data/{0}-cal-sec-star/hists-{1}.npy'.format(name, suffix))
			for i in range(detector_size[0]):
				for j in range(detector_size[1]):
					hist_list[i*detector_size[1]+j] += hists[i*detector_size[1]+j]
		
		centroids = []
		for i in range(detector_size[0]):
			for j in range(detector_size[1]):
				H = hist_list[i*detector_size[1]+j]    
				t = 1
				print 'threshold:{0}'.format(t) 
				#plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'))
				#plt.show()
				data = H.byteswap(True).newbyteorder()
				data = data.copy(order='C')
				data = data.byteswap(True).newbyteorder()
				c_data = c3.find_centroid(data, t)
				if c_data is not None:
					cy, cx, max_value, flux = c_data
					centroids.append([-0.2+cx*0.4/hist_size[0], -0.2+cy*0.4/hist_size[0]])
				else:
					centroids.append([0, 0])
		name = '{0}-{1}-{2}'.format(int(name_list[0].split('_')[3]), int(name_list[-1].split('_')[3]), suffix)

		np.save('../data/distortion/{0}-centroids.npy'.format(name), centroids)

		c = np.load('../data/distortion/{0}-centroids.npy'.format(name))
		cen = c.reshape((50,50,2))
		cx = np.zeros((50,50))
		cy = np.zeros((50,50))
		for i in range(50):
			for j in range(50):
				cx[i,j] = cen[i,j,0]
				cy[i,j] = cen[i,j,1]
		plt.imshow(cx,interpolation='None', cmap=plt.get_cmap('Greys'), vmin=-0.05, vmax=0.05)
		#plt.title(r'{0} $\Delta x$'.format(name))
		plt.title(r'$\Delta x\  median\ XA=25$')
		cbar = plt.colorbar()
		cbar.set_label(r'$\Delta x$ [pixel]', rotation=270)
		plt.tight_layout()
		#plt.savefig('/home/dw1519/galex/plots/cx/{0}-cx_new.png'.format(name))
		plt.savefig('/home/dw1519/galex/plots/{0}-cx_new.png'.format(name), dpi=190)
		plt.clf()
		plt.imshow(cy,interpolation='None', cmap=plt.get_cmap('Greys'), vmin=-0.05, vmax=0.05)
		#plt.title(r'{0} $\Delta y$'.format(name))
		plt.title(r'$\Delta y\ median\ XA=25$')
		cbar = plt.colorbar()
		cbar.set_label(r'$\Delta y$ [pixel]', rotation=270)
		plt.tight_layout()
		#plt.savefig('/home/dw1519/galex/plots/cy/{0}-cy_new.png'.format(name))
		plt.savefig('/home/dw1519/galex/plots/{0}-cy_new.png'.format(name), dpi=190)
		plt.clf()