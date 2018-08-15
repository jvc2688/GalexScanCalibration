import numpy as np
from astropy.io import fits as pyfits
from numpy.random import multivariate_normal
from numpy.random import uniform
import imagetools
import c3
import matplotlib.pyplot as plt


def hist_g0(data, wcs, imsz):
  foc = wcs.all_world2pix(data,1)
  H,xedges,yedges=np.histogram2d(foc[:,1], foc[:,0],\
                            bins=imsz, range=[ [0.5, imsz[0]+0.5],[0.5,imsz[1]+0.5] ])
  return H

def hist_g(data, wcs, imsz):
  foc = wcs.all_world2pix(data,1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

def hist(data, pixsz, skyrange):
    imsz = np.array(skyrange)/pixsz
    H,xedges,yedges=np.histogram2d(data[:,0], data[:,1],\
                            bins=imsz, range=([ [-skyrange[0]/2.,skyrange[0]/2.],[-skyrange[1]/2.,skyrange[1]/2.] ]))
    return H

if __name__ == '__main__':

    if True:
        x = uniform(-0.495,0.495, 200).astype(float)
        y = uniform(-0.495,0.495, 200).astype(float)
        pos_list = np.array([x,y]).T
        np.savetxt('../fits/scan_map/source_new.dat', pos_list, fmt='%.8f')
        #pos_list = [[0,0],[0.1, 0.1],[0.11, -0.1],[-0.05,-0.01],[0.08,0.08],[-0.12,0.02]]
        fwhm = (3./3600.)**2
        cov = [[fwhm, 0], [0, fwhm]]

        skypos = [0,0]
        skyrange = [1., 1.]
        pixsz = 0.000416666666666667
        imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz).astype(int)
        print imsz
        wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
        count = np.zeros(imsz)
        photon_list = []
        offsets = []
        i=0
        for pos in pos_list:
            print pos
            data= multivariate_normal(pos, cov, 100000)
            print np.mean(data, axis=0)
            photon_list.append(data)
            off = data-pos
            offsets.append(data-pos)
            print data.shape
            h = hist_g(data, wcs, imsz)
            count += h
            fig = plt.gcf()
            fig.set_size_inches(8,8)
            plt.axhline(y=0, color='r')
            plt.axvline(x=0, color='r')
            plt.plot(off[:,0]*3600, off[:,1]*3600, 'o', alpha=0.01)
            plt.xlim(-0.003*3600, 0.003*3600)
            plt.ylim(-0.003*3600, 0.003*3600)
            plt.savefig('../fits/scan_map/plot/{0}.png'.format(i))
            plt.clf()
            i+=1
        offset = np.concatenate(offsets, axis=0)
        np.save('../fits/scan_map/photon_list.npy', photon_list)
        tranges = [[0,0]]
        hdu = pyfits.PrimaryHDU(count)
        hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('../fits/scan_map/fake_new_small.fits', clobber=False)

    if False:
        x = [0.0125]
        y = [-0.00812]
        pos_list = np.array([x,y]).T
        fwhm = (5./3600.)**2
        cov = [[fwhm, 0], [0, fwhm]]

        skypos = [0,0]
        skyrange = [0.05, 0.05]
        pixsz = 0.000416666666666667
        pixsz = 0.0002
        imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz).astype(int)
        print imsz
        wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
        count = np.zeros(imsz)
        photon_data = None
        for pos in pos_list:
            photon_data = multivariate_normal(pos, cov, 50000)
            print photon_data.shape
            h = hist_g(photon_data, wcs, imsz)
            count += h
        '''
        tranges = [[0,0]]
        hdu = pyfits.PrimaryHDU(count)
        hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('../fits/scan_map/fake2.fits', clobber=False)
        '''
        t = 1

        H = h.copy()
        data = H.byteswap(True).newbyteorder()
        data = data.copy(order='C')
        data = data.byteswap(True).newbyteorder()
        c_data = c3.find_centroid(data, t)
        if c_data is not None:
            cy, cx, max_value, flux = c_data
            cy+=0#0.5
            cx+=0#0.5
            centroid = wcs.wcs_pix2world(wcs.sip_foc2pix([[cx, cy]],1),1)[0]
        else:
            centroid = [0.,0.]
            max_value = 0
            flux = 500
        if centroid[0]>1:
            centroid[0] = centroid[0]-360.
        #print 'max:{0}'.format(max_value)
        fig = plt.gcf()
        fig.set_size_inches(8,8)
        plt.imshow(h,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[0.025,-0.025,0.025,-0.025], origin='upper', imlim=[-0.02,0.02, -0.02,0.02])
        #plt.imshow(hp,interpolation='None', cmap=plt.get_cmap('Greys'), extent=[0.02,-0.02,0.02,-0.02])
        plt.xlabel(r'$\Delta RA$')
        plt.ylabel(r'$\Delta Dec$')
        plt.plot(centroid[0],centroid[1],'+r', markersize=8)
        plt.xlim(0.02,-0.02)
        plt.ylim(0.02,-0.02)
        plt.tight_layout()
        plt.show()
        #plt.subplots_adjust(left=0.12, bottom=0.07, right=0.95, top=0.97, wspace=0, hspace=0)
        #plt.show()
        plt.clf()

        hp = hist(photon_data, pixsz, skyrange)
        plt.imshow(hp,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[0.025,-0.025,0.025,-0.025], origin='upper', imlim=[-0.02,0.02, -0.02,0.02])
        plt.show()
        data = hp.byteswap(True).newbyteorder()
        data = data.copy(order='C')
        data = data.byteswap(True).newbyteorder()
        c_data = c3.find_centroid(data, t)
        cy, cx, max_value, flux = c_data
        print -skyrange[0]/2.+cy*pixsz, -skyrange[1]/2.+cx*pixsz
        print centroid, max_value, flux
