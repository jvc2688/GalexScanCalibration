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
        num = int(re.split('_', name)[3])
        scan = '%04d'%num
        print scan
        catalog = '/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
        if os.path.isfile(catalog):
        	
