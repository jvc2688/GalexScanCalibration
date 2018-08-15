import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import imagetools
from matplotlib.colors import LogNorm
import re
import glob
import os
import gc
import sys

def get_key(path):
	return int(re.split('/|\.', path)[-2])

if __name__ == '__main__':
	if True:
		sub = sys.argv[1]
		with open(sub,'r') as f:
			name_list = f.read().splitlines()
		for name in name_list:
			#scan = re.split('_', name)[3]
			qfile = '/scratch/dw1519/galex/co/remove/'+name+'.npy'
			if os.path.isfile(qfile):
				qhists = np.load(qfile)
				print 'existed'
			else:
				photon_list =  sorted(glob.glob("/beegfs/dw1519/galex/data/"+name+"/split/*.csv"), key=get_key)
				#print photon_list
				qhists = []
				data_list = []
				#bins = np.arange(33)
				for f in photon_list:
					print f
					data = np.loadtxt(f,delimiter=',', usecols=[0,3,4,5,6,7])
					data_list.append(data)
					if len(data.shape)<2:
						continue
					mask = (data[:,2]>=2) & (data[:,3]>5)
					qhists.append([data[0,0], float(np.sum(mask))/data.shape[0]])
				np.save(qfile, qhists)
				data = np.concatenate(data_list,axis=0)
				np.save('/scratch/dw1519/galex/data/photons/'+name+'.npy', data)
				data=data_list=None
				qhists = np.array(qhists)
			print qhists.shape
			plt.plot((qhists[:,0]-qhists[0,0])/1000., qhists[:,1])
			plt.title('{0}'.format(name)) 
			plt.xlabel(r't [s]')
			#plt.savefig(out+'/xya_{0}.png'.format((i+1)*100), dpi=190)
			plt.savefig('/scratch/dw1519/galex/co/remove/'+name+'.pdf', dpi=190)
			plt.clf()
			qhists=None
			gc.collect()
