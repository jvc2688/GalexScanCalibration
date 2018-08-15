import numpy as np 
from astropy.io import fits as pyfits
import re
import glob
from astropy.nddata.utils import block_reduce
from astropy.wcs import WCS


def get_resampled_wcs(wcs, factor, downsampled):
    """
    Get resampled WCS object.
    """
    wcs = wcs.deepcopy()

    if not downsampled:
        factor = 1. / factor

    wcs.wcs.cdelt *= factor
    wcs.wcs.crpix = (wcs.wcs.crpix - 0.5) / factor + 0.5
    return wcs


def in_map(path, name):
    count = pyfits.open(path+'/count_map_%s_count.fits'%name)
    exp = pyfits.open(path+'/count_map_%s_exp.fits'%name)

    cm = count[0].data
    em = exp[0].data

    em[np.where(em<0.00001)]=1.
    print np.min(em)
    print np.max(em)

    ratio = cm/em
    print np.min(cm)
    print np.max(cm)

    print np.min(ratio)
    print np.max(ratio)

    hdu = pyfits.PrimaryHDU(ratio)
    hdu.header = count[0].header
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(path+'/count_map_%s_in.fits'%name, clobber=False)



if __name__ == '__main__':

    if True:
        name_list = glob.glob("../name_new/p2/*")
        name_list = ['name2057-2102']
        name_list = ['name_185', 'name_239', 'name_248', 'name_221','name_230', 'name_3029', 'name_3038', 'name_3047', 'name_3245', 'name_3299', 'name_3317', 'name_3488']
        name_list = ['3209', '0014']
        name_list = ['0023', '0032', '0050', '0356', '0392', '0743', '1103', '2381', '3587']
        name_list = ['5-23']
        name_list = ['0014', '0032', '0059', '0203', '0239']
        name_list = ['0239']
        name_list = ['0023', '0032', '0203', '0239', '0446', '0464', '0473', '0806', '0815', '1301', '1310', '1319',\
                    '1616', '1634', '1679', '2174', '2183', '2192', '2714', '2750', '3236', '3245', '3281']
        name_list = ['2561']
        name_list = ['0050']
        date='11-07-2017'
        date1='08-21-2017/0001-0500'
        factor = 4
        for name in name_list:
            #name = re.split('/', name)[3]
            if name == 'name2444-2480':
                continue
            print name
            count = pyfits.open('../fits/scan_map/'+date1+'/count_map_%s_count.fits'%name)
            exp = pyfits.open('../fits/scan_map/'+date+'/count_map_%s_exp.fits'%name)

            cm = block_reduce(count[0].data, factor)
            em = block_reduce(exp[0].data, factor, np.mean)

            w = WCS('../fits/scan_map/'+date1+'/count_map_%s_count.fits'%name)
            wn = get_resampled_wcs(w,factor,True)
            header = wn.to_header()

            hdu = pyfits.PrimaryHDU(em)
            hdu.header = count[0].header
            hdu.header['NAXIS1'] = em.shape[1]#header['NAXIS1']
            hdu.header['NAXIS2'] = em.shape[0]#header['NAXIS2']
            hdu.header['CRPIX1'] = header['CRPIX1']
            hdu.header['CRPIX2'] = header['CRPIX2']
            hdu.header['CDELT1'] = header['CDELT1']
            hdu.header['CDELT2'] = header['CDELT2']
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto('../fits/scan_map/'+date+'/count_map_%s_ds_exp.fits'%name, clobber=False)

            em[np.where(em<0.001)]=1.
            print np.min(em)
            print np.max(em)

            ratio = cm/em
            print np.min(cm)
            print np.max(cm)

            print np.min(ratio)
            print np.max(ratio)

            w = WCS('../fits/scan_map/'+date1+'/count_map_%s_count.fits'%name)
            wn = get_resampled_wcs(w,factor,True)
            header = wn.to_header()

            hdu = pyfits.PrimaryHDU(ratio)
            hdu.header = count[0].header
            hdu.header['NAXIS1'] = ratio.shape[1]#header['NAXIS1']
            hdu.header['NAXIS2'] = ratio.shape[0]#header['NAXIS2']
            hdu.header['CRPIX1'] = header['CRPIX1']
            hdu.header['CRPIX2'] = header['CRPIX2']
            hdu.header['CDELT1'] = header['CDELT1']
            hdu.header['CDELT2'] = header['CDELT2']
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto('../fits/scan_map/'+date+'/count_map_%s_ds_in.fits'%name, clobber=False)

    if False:
        count = pyfits.open('../fits/count_map_name671_gPr_cata_10_gal.fits')
        exp = pyfits.open('../fits/exp/count_map_name671_gPr_flip_gal_dead_sec.fits')

        cm = count[0].data
        em = exp[0].data

        em[np.where(em<0.00001)]=1.
        print np.min(em)
        print np.max(em)

        ratio = cm/em
        print np.min(cm)
        print np.max(cm)

        print np.min(ratio)
        print np.max(ratio)

        hdu = pyfits.PrimaryHDU(ratio)
        hdu.header = count[0].header
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('../fits/count_map_name671_gPr_cata_10_corr_gal.fits', clobber=False)

    #focal plane
    if False:
        focal = pyfits.open('count_map_05-68_gPr_cata_10_device_r0_nostar.fits')
        flat = pyfits.open('NUV_flat.fits')

        focal_m = focal[0].data
        flat_m = flat[0].data

        #flat_m = np.flipud(np.rot90(flat_m))
        
        flat_m = np.swapaxes(flat_m,0,1)

        print np.median(focal_m)*0.85, np.median(focal_m)*1.15
        print np.median(flat_m)*0.85, np.median(flat_m)*1.15

        #flat_m[np.where(flat_m<0.00001)]=1.
        flat_m[np.where(flat_m<0.001)]=1.

        #print np.min(flat_m)
        #print np.max(flat_m)

        ratio = focal_m/flat_m
        #print np.min(focal_m)
        #print np.max(focal_m)

        #print np.min(ratio)
        #print np.max(ratio)

        print np.median(ratio)*0.85, np.median(ratio)*1.15


        hdu = pyfits.PrimaryHDU(ratio)
        hdu.header = focal[0].header
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('count_map_05_68_gPr_cata_10_device_r0_ratio_hi.fits', clobber=False)

