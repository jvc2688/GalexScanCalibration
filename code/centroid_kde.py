import numpy as np 
from astropy.io import fits as pyfits
from sklearn.neighbors.kde import KernelDensity
import csv
import math

def load_data(filename):
	data = []
	with open(filename, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			if math.fabs(float(row[0]))>0.008 or math.fabs(float(row[1]))>0.008:
				pass
			else:
				data.append(row)
		csvfile.close()
	return data

def find_centroid(data, bandwidth=0.003, iter_num=6, halfwidth=0.02):
	kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(data)
	grid = 10
	position = np.array([0,0])
	#halfwidth = 0.02
	for i in range(iter_num):
		low = position-halfwidth
		high = position+halfwidth
		X, Y = np.mgrid[low[0]:high[0]:10j, low[1]:high[1]:10j]
		positions = np.vstack([X.ravel(), Y.ravel()]).T
		img = kde.score_samples(positions)
		position = positions[np.argmax(img)]
		halfwidth = halfwidth*2./(grid-1.)
	return position

if __name__ == '__main__':
	pass


