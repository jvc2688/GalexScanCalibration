import centroid_plot_csv_new_0 as centroid_csv
import numpy as np
import sys
import os

if __name__ == '__main__':
	sub_name = sys.argv[1]
	print sub_name
	scan_list = []
	with open(sub_name, 'r') as f:
		scan_list = f.read().splitlines()
	for name in scan_list:
		print name
		asprta = '/scratch/dw1519/galex/AIS_GAL_SCAN/asprta/'+name+'-asprta.fits'
		out_suffix = 'cal-sec'
		tmp_dir = '/scratch/dw1519/galex/data/'+name+'/cata'
		p_num=0
		if os.path.isfile('/scratch/dw1519/galex/plots/0/%s-%s_ra.pdf'%(name, out_suffix)):
			print 'skip'
		else:
			centroid_csv.generate_new_offsets_new(name, asprta, out_suffix, tmp_dir, p_num)