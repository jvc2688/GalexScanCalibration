import sys
import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits as pyfits
import centroid_plot

if __name__ == '__main__':
	initial = 300
	final = 1342
	length = final-initial+1
	offsets = centroid_plot.get_centers(initial, final)

	for i in range(length-1):
		offsets[i+1] = offsets[i]+offsets[i+1]

	trange = np.arange(final-initial+1)+initial
	'''
	plt.plot(trange, offsets[:,0],'-b')
	plt.xlabel('time/s')
	plt.ylabel('RA/degree')
	plt.savefig('../plots/ra%d_%d_cul.png'%(initial,final),dpi=190)
	plt.clf()
	plt.plot(trange, offsets[:,1],'-b')
	plt.xlabel('time/s')
	plt.ylabel('DEC/degree')
	plt.savefig('../plots/dec%d_%d_cul.png'%(initial,final),dpi=190)
	'''
	filename= '../data/cal/AIS_GAL_SCAN_00005_0001-asprta_cul_b.fits'
	asprate_file  = pyfits.open(filename, mode='update')
	data = asprate_file[1].data
	data['ra'][initial+1:initial+1+length] = data['ra'][initial+1:initial+1+length] - offsets[:,0]
	data['dec'][initial+1:initial+1+length] = data['dec'][initial+1:initial+1+length] - offsets[:,1]
	asprate_file.flush()
	asprate_file.close()
