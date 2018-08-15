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
	if True:
		suffix_list = ['q_0', 'q_1', 'q_2']
		#suffix_list = ['xa_0', 'xa_1', 'xa_2']

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
			suffix = suffix_list[0]
			name = '{0}-{1}-{2}'.format(int(name_list[0].split('_')[3]), int(name_list[-1].split('_')[3]), suffix)
			cen0 = np.load('../data/distortion/{0}-centroids.npy'.format(name)).reshape((50,50,2))
			for suffix in suffix_list:	
				name = '{0}-{1}-{2}'.format(int(name_list[0].split('_')[3]), int(name_list[-1].split('_')[3]), suffix)
				c = np.load('../data/distortion/{0}-centroids.npy'.format(name))
				cen = c.reshape((50,50,2))
				cx = np.zeros((50,50))
				cy = np.zeros((50,50))
				for i in range(50):
					for j in range(50):
						cx[i,j] = cen[i,j,0]
						cy[i,j] = cen[i,j,1]
				plt.imshow(cx-cen0[:,:,0],interpolation='None', cmap=plt.get_cmap('Greys'), vmin=-0.05, vmax=0.05)
				plt.title(r'$\Delta x$ difference'.format(name))
				cbar = plt.colorbar()
				cbar.set_label(r'$\Delta x$', rotation=270)
				plt.tight_layout()
				plt.savefig('/home/dw1519/galex/plots/dif1_{0}-cx_new.png'.format(name), dpi=190)
				plt.clf()
				plt.imshow(cy-cen0[:,:,1],interpolation='None', cmap=plt.get_cmap('Greys'), vmin=-0.05, vmax=0.05)
				plt.title(r'$\Delta y$ difference'.format(name))
				cbar = plt.colorbar()
				cbar.set_label(r'$\Delta y$', rotation=270)
				plt.tight_layout()
				plt.savefig('/home/dw1519/galex/plots/dif1_{0}-cy_new.png'.format(name), dpi=190)
				plt.clf()
				cen0 = cen

	if False:
		suffix_list = ['q_low', 'q_high']
		scan_list = ['5-185', '194-437', '446-887', '896-1157', '1166-1373']

		for suffix in suffix_list:
			#suffix = suffix_list[0]
			name = '{0}-{1}'.format(scan_list[0], suffix)
			cen0 = np.load('../data/distortion/{0}-centroids.npy'.format(name)).reshape((50,50,2))
			for name in scan_list:	
				name = '{0}-{1}'.format(name, suffix)
				c = np.load('../data/distortion/{0}-centroids.npy'.format(name))
				cen = c.reshape((50,50,2))
				cx = np.zeros((50,50))
				cy = np.zeros((50,50))
				for i in range(50):
					for j in range(50):
						cx[i,j] = cen[i,j,0]
						cy[i,j] = cen[i,j,1]
				plt.imshow(cx-cen0[:,:,0],interpolation='None', cmap=plt.get_cmap('Greys'), vmin=-0.05, vmax=0.05)
				plt.title(r'{0} $\Delta x$ difference'.format(name))
				cbar = plt.colorbar()
				cbar.set_label(r'$\Delta x$', rotation=270)
				plt.savefig('/home/dw1519/galex/plots/cx/dift_{0}-cx_new.png'.format(name), dpi=190)
				plt.clf()
				plt.imshow(cy-cen0[:,:,1],interpolation='None', cmap=plt.get_cmap('Greys'), vmin=-0.05, vmax=0.05)
				plt.title(r'{0} $\Delta y$ difference'.format(name))
				cbar = plt.colorbar()
				cbar.set_label(r'$\Delta y$', rotation=270)
				plt.savefig('/home/dw1519/galex/plots/cy/dift_{0}-cy_new.png'.format(name), dpi=190)
				plt.clf()
				cen0 = cen
