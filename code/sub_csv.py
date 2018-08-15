import sys
from subprocess import call

if __name__ == '__main__':
	sub_name = sys.argv[1]
	print sub_name
	scan_list = []
	with open(sub_name, 'r') as f:
		scan_list = f.read().splitlines()
	for name in scan_list:
		print name
		call(["/home/dw1519/dw1519/anaconda/bin/python", "csv2bin.py", "/scratch/dw1519/galex/data/", name])
