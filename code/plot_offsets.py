import matplotlib.pyplot as plt
import numpy as np 
import sys

if __name__ == '__main__':
	name = sys.argv[1]
	offsets = np.load('../data/%s/cata/offsets_inter_half_sec.npy'%name)
	plt.plot(offsets[:,0], '.b', markersize=0.5)
	plt.show()