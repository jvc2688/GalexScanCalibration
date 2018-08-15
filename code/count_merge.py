import numpy as np 
from astropy.io import fits as pyfits
import re
import glob
import gc
import sys
from itertools import groupby


def get_key(name):
    num = int(re.split('_', name)[3])
    return '%04d'%num



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
        name_list = ['0284', '0293', '0302', '0311', '0320', '0329', '0338', '0347', '0356']
        name_list = ['3596', '0005', '0014', '0023', '0032', '0041', '0050', '0059', '0068', '0086', '0095', '0104']
        name_list = ['1724', '1733', '1742','1751', '1760', '1769', '1778', '1787', '1796', '1805']
        name_list = ['0185', '0194', '0203','0212', '0221', '0230', '0239', '0248', '0257', '0284']
        name_list = ['0086', '0095', '0104', '0113', '0122', '0140', '0149', '0158', '0167', '0176', '0185', '0194']
        name_list = ['2345', '2354', '2363', '2372', '2381', '2390', '2399', '2408', '2417', '2426', '2435', '2444']

        sky = sys.argv[1]
        print sky
        with open('/scratch/dw1519/galex/count_pbs/'+sky, 'r') as f:
            scan_list = f.read().splitlines()

        name_list = []
        for key, group in groupby(scan_list, get_key):   
            name_list.append(key)
        print name_list

        #sky = '9-19'
        i=0
        for name in name_list:
            #name = re.split('/', name)[3]
            if name == 'name2444-2480':
                continue
            print name
            fod = '270-280'
            out = 'fix'
            if i==0:
                count = pyfits.open('../fits/scan_map/12-05-2017/'+fod+'/count_map_%s_count.fits'%name)
                exp = pyfits.open('../fits/scan_map/12-05-2017/'+fod+'/count_map_%s_exp.fits'%name)
                c = count[0].data
                e = exp[0].data
                c[np.isnan(c)] = 0
                c[c>5000] = 0
                c[c<0] = 0
                e[np.isnan(e)] = 0
                e[e>5000] = 0
                e[e<0] = 0
                cm = c
                em = e
                print np.max(cm)
                print np.max(em)
            else:
                ch = pyfits.open('../fits/scan_map/12-05-2017/'+fod+'/count_map_%s_count.fits'%name)
                eh = pyfits.open('../fits/scan_map/12-05-2017/'+fod+'/count_map_%s_exp.fits'%name)
                c = ch[0].data
                e = eh[0].data
                c[np.isnan(c)] = 0
                c[c>5000] = 0
                c[c<0] = 0
                e[np.isnan(e)] = 0
                e[e>5000] = 0
                e[e<0] = 0
                cm += c
                em += e
                ch.close()
                eh.close()
                print np.max(cm)
                print np.max(em)
                c=None
                e=None
            i+=1
            gc.collect()


        hdu = pyfits.PrimaryHDU(cm.astype('float32'))
        hdu.header = count[0].header
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('../fits/scan_map/12-05-2017/'+out+'/count_map_%s_count.fits'%sky, clobber=False)


        hdu = pyfits.PrimaryHDU(em.astype('float32'))
        hdu.header = exp[0].header
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('../fits/scan_map/12-05-2017/'+out+'/count_map_%s_exp.fits'%sky, clobber=False)
        '''

        count = pyfits.open('../fits/scan_map/12-05-2017/count_map_%s_count.fits'%sky)
        cm = count[0].data

        exp = pyfits.open('../fits/scan_map/12-05-2017/count_map_%s_exp.fits'%sky)
        em = exp[0].data
        '''
        em[np.where(em<0.00001)]=1.
        print np.min(em)
        print np.max(em)

        ratio = cm/em
        print np.min(cm)
        print np.max(cm)

        print np.min(ratio)
        print np.max(ratio)

        hdu = pyfits.PrimaryHDU(ratio.astype('float32'))
        hdu.header = count[0].header
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('../fits/scan_map/12-05-2017/'+out+'/count_map_%s_in.fits'%sky, clobber=False)

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

