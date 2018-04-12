from astropy.io import fits as pyfits
import numpy as np

def main():
	for i in range(0,41):
		name_file = 'name%d'%i
		with open('../name/%s'%name_file) as f:
			name_list = f.read().splitlines()
		#print name_list
		for name in name_list:
			hdu = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
			data = hdu[1].data
			roll = np.median(data['Roll'])
			if roll<330 and roll >315:
				print '%s %f'%(name, roll)

if __name__ == '__main__':
	main()