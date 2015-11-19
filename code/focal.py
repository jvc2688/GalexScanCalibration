import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits as pyfits

hdulist = pyfits.open('../../gPhoton/Galex/AIS_GAL_SCAN_00005_0001-nd-raw6.fits')

x = hdulist[1].data['phb1']
y = hdulist[1].data['phb2']

xedges = np.arange(256)
yedges = np.arange(256)
print x.dtype, x.shape
print y

H, xedges, yedges = np.histogram2d(y, x, bins=(xedges, yedges))

print H
print H.shape

fig = plt.figure()
axes = fig.add_subplot(111)
im = axes.imshow(H, cmap='gray_r', aspect='auto', interpolation="none")

plt.savefig('focal.png', dpi=190)
plt.clf()



