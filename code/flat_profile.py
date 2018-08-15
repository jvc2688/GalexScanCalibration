import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
from astropy.io import fits as pyfits

def convert1(flat):
    r=[]
    mask = np.isfinite(flat)
    for i in range(800):
        rs = flat[:,i][mask[:,i]]
        if len(rs)>0:
            r.append(np.mean(rs))
        else:
            r.append(0)
    return np.array(r)


if __name__ == '__main__':
    bins = [0.,0.015,0.02,0.023,0.026,0.03,0.1]
    flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
    flat0 = flat_hdu[0].data
    profile0 = convert1(flat0)
    mask = profile0>0
    #profile0[~mask] = 1.
    ratio = profile0/profile0
    plt.plot(profile0,'-', label='pipeline')
    for i in range(6):        
        flat = np.load('/scratch/dw1519/galex/data/star_photon/super_iterate{0}/flat.npy'.format(i+4))
        profile = convert1(flat)/0.6*0.75
        mask2 = profile>0
        #profile[~mask2] = 1.
        plt.plot(profile,'-',label='{0}-{1}'.format(bins[i],bins[i+1]))
    plt.legend()
    #plt.ylim(0.5,1.)
    plt.xlabel('x [pixel]')
    plt.ylabel('Average sensitivity')
    plt.savefig('/scratch/dw1519/galex/data/star_photon/profile_x.png', dpi=190)
    plt.clf()

    flat0 = rotate(flat0,-90.,reshape=False,order=1,prefilter=False)
    profile0 = convert1(flat0)
    mask = profile0>0
    #profile0[~mask] = 1.
    ratio = profile0/profile0
    plt.plot(profile0,'-', label='pipeline')

    for i in range(6):        
        flat = np.load('/scratch/dw1519/galex/data/star_photon/super_iterate{0}/flat.npy'.format(i+4))
        flat = rotate(flat,-90.,reshape=False,order=1,prefilter=False)
        profile = convert1(flat)/0.6*0.72
        mask2 = profile>0
        #profile[~mask2] = 1.
        plt.plot(profile,'-',label='{0}-{1}'.format(bins[i],bins[i+1]))
    plt.legend()
    #plt.ylim(0.5,1.)
    plt.xlabel('y [pixel]')
    plt.ylabel('Average sensitivity')
    plt.savefig('/scratch/dw1519/galex/data/star_photon/profile_y.png', dpi=190)
    plt.clf()