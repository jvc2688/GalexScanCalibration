import csv
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import re
import os
from itertools import groupby
from operator import itemgetter

def get_key(path):
	name = os.path.basename(path)
	print re.split('-|\.', name)
	return int(re.split('-|\.', name)[1])

def get_group_key(arg):
	name = os.path.basename(arg[1])
	#print re.split('-|\.', name)
	index = arg[0]
	return index-int(re.split('-|\.', name)[1])

if __name__ == '__main__':
	name = 'AIS_GAL_SCAN_00032_0001'

	output = '../data/%s/star/extra/bkg/star.csv'%(name)
	print output
	path = os.path.dirname(output)
	if not os.path.exists(path):
		os.makedirs(path)
	else:
		print 'exists'
	for i in range(5989):
		star_num = '%d'%i#sys.argv[1]
		csv_list = glob.glob("../data/%s/star/extra/star_bkg_%s-*.csv"%(name, star_num))
		csv_list = sorted(csv_list, key=get_key)
		print csv_list
		groups = []
		for k, g in groupby(enumerate(csv_list), get_group_key):
			groups.append(map(itemgetter(1), g))
		print groups
		plot_num = 0
		for csvs in groups:
			q_count = []
			count = []
			time_list = []
			fout = open('../data/%s/star/extra/bkg/star_bkg_%s-%d.csv'%(name, star_num, plot_num), 'ab')
			outwriter = csv.writer(fout)
			last = 0
			for csv_file in csvs:
				with open(csv_file, 'rb') as file:
					reader = csv.reader(file)
					#last = int(float(reader.next()[0])/1000.)
					#time_list.append(last)
					#count.append(1)
					for row in reader:
						outwriter.writerow(row)
						'''
						time = float(row[0])#time = int(float(row[0])/1000.)
						if time - last > 1000:
							count.append(1)
							last = time
							time_list.append(last/1000.)
						'''
						'''
							count.append(1)
							q_count.append(float(row[5]))
							last=time
							time_list.append(last)
						'''
						'''
						else:
							count[-1] += 1
							#count[-1] += 1
							#q_count[-1] += float(row[5])
						'''
			#q = np.array(q_count)/np.array(count)
			fout.close()
			'''
			plt.plot(time_list, count)
			plt.xlabel('t/s')
			plt.ylabel('Number of Photons')
			plt.savefig('../plots/%s/star/extra/%s-%d.png'%(name, star_num, plot_num), dpi=190)
			plt.clf()
			'''
			plot_num +=1
