import numpy as np
import sys
import os
import centroid_plot_csv_new_test as centroid_csv
import re

if __name__ == '__main__':
    sub_name = sys.argv[1]
    print sub_name
    scan_list = []
    with open(sub_name, 'r') as f:
        scan_list = f.read().splitlines()
    for name in scan_list:
        print name
        num = int(re.split('_', name)[3])
        scan = '%04d'%num
        print scan
        asprta = '/scratch/dw1519/galex/AIS_GAL_SCAN/asprta/'+name+'-asprta.fits'
        tmpdir='/scratch/dw1519/galex/data/'+name+'/cata'
        suffix='cal-sec'
        centroid_csv.generate_new_offsets_new(name, asprta, suffix, tmpdir, 10)
