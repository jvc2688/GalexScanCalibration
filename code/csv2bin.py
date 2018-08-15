import numpy as np
import sys
import glob
import re
import os
import subprocess

def get_key(path):
	return int(re.split('/|\.',path)[7])

if __name__ == '__main__':
	if False:
		data_path = sys.argv[1]
		name = sys.argv[2]
		if len(sys.argv)>3:
			suffix = sys.argv[3]
		else:
			suffix = ''
		out_path = '/scratch/dw1519/galex/data/'
		csvs = sorted(glob.glob(data_path+'{0}{1}/split/*.csv'.format(name,suffix)), key=get_key)
		data_list=[]#np.array([[],[],[],[],[],[]]).T
		for csv in csvs:
			print csv
			data = np.loadtxt(csv,delimiter=',', usecols=[0,3,4,5,6,7])
			if len(data.shape)<2:
				print 'skip'
				continue
			data_list.append(data)
		data = np.concatenate(data_list,axis=0)
		np.save(out_path+'{0}/photon_list.npy'.format(name), data)

	if True:
		data_path = sys.argv[1]
		name = sys.argv[2]
		if len(sys.argv)>3:
			suffix = sys.argv[3]
		else:
			suffix = ''
		print suffix
		out_path = '/beegfs/dw1519/galex/data/'
		output = out_path+'{0}{1}/test'.format(name,suffix)
		dir = os.path.dirname(output)
		if not os.path.exists(dir):
			os.makedirs(dir)
		if not os.path.isfile(out_path+'{0}{1}/split/lookup.pkl'.format(name, suffix)):
			subprocess.call(["cp", '-r', data_path+'{0}{1}/split'.format(name,suffix), out_path+'{0}{1}/.'.format(name, suffix)])
		else:
			print 'skip'

