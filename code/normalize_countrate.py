from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import load_data
from matplotlib.colors import LogNorm
from astropy import wcs
import glob
import os
import re
from scipy.ndimage.interpolation import rotate
#plt.style.use('seaborn-paper')

def convert1(flat):
    r=[]
    mask = np.isfinite(flat)&(flat>0.)
    for i in range(800):
        rs = flat[:,i][mask[:,i]]
        if len(rs)>0:
            r.append(np.mean(rs))
        else:
            r.append(0)
    return np.array(r)

def plot(date, scan):
    #date = '08-21-2017'
    print(scan)
    tycho2 = pyfits.open('/scratch/dw1519/galex/data/tycho2.fits')[1].data
    cfile='/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
    df = load_data.load_catalog(cfile)
    c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
    '''
    if scan=='3488':
        cut = df['gb']>1.
        c = c[cut]
        df = df[cut]
    '''
    print(c.shape)
    catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
    mask=d2d<0.001*u.degree
    N = np.sum(mask)
    print(N)

    f = plt.figure(figsize=(12,10))
    ax = plt.subplot(2, 2, 1)
    ax.grid(False)
    plt.scatter((c[mask].l.deg-catalog[idx[mask]].l.deg)*3600,\
        (c[mask].b.deg-catalog[idx[mask]].b.deg)*3600,\
         c=df['nuv'][mask], alpha=0.2, edgecolors='face', cmap=plt.get_cmap('jet'))
    plt.axis('equal')
    plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
    plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
    #plt.title('0023 NUV<15, star={0} {1:.0f}%, photon={2} {3:.0f}%'.format(len(bright), len(bright)/total_star*100,\
    #                                                len(bright_co), len(bright_co)/total_photon*100))
    plt.xlim(-3.0,3.0)
    plt.ylim(-3.0,3.0)
    plt.ylabel('gb(SEx-T2)')
    plt.title('gl(SEx-T2), N={0}'.format(N))
    plt.colorbar()


    ais = pyfits.open('/scratch/dw1519/galex/data/MyTable_4_jvc2688.fit')[1].data
    ais_cat = SkyCoord(ais['glon']*u.degree, ais['glat']*u.degree, frame='galactic')
    idx, d2d, d3d = c.match_to_catalog_sky(ais_cat)  
    mask=d2d<0.001*u.degree
    N = np.sum(mask)
    print('ais match:{0}'.format(N))
    plt.subplot(2, 2, 2)
    plt.scatter(ais['nuv'][idx[mask]],\
        df['nuv'][mask]-ais['nuv'][idx[mask]],\
         c=ais['glat'][idx[mask]], alpha=0.5, edgecolors='face', s=0.2, cmap=plt.get_cmap('jet'), vmin=-10, vmax=10)
    plt.axhline(y=0, color='k', linestyle='dashed', linewidth=1)
    nuv_list = [[14,15],[15,16],[16,17],[17,18],[18,19],[19,20]]
    mask_list = [(e[0]<=ais['nuv'][idx[mask]]) & (e[1]>ais['nuv'][idx[mask]]) for e in nuv_list]
    x = [np.median(ais['nuv'][idx[mask]][m]) for m in mask_list]
    y = [np.median((df['nuv'][mask]-ais['nuv'][idx[mask]])[m]) for m in mask_list]
    ye = [np.std((df['nuv'][mask]-ais['nuv'][idx[mask]])[m]) for m in mask_list]
    plt.errorbar(x, y, yerr=ye)
    plt.xlim(12,22)
    plt.ylabel('NUV (SEx-AIS)')
    plt.title('NUV AIS, N={0}'.format(N))
    plt.ylim(-1,1)
    plt.colorbar()

    plt.subplot(2, 2, 3)
    plt.hist(df['nuv'],bins=18,range=(12,22))
    plt.title('NUV SExtractor')
    plt.xlim(12,22)

    plt.subplot(2, 2, 4)
    plt.plot(df['nuv'], df['FWHM_IMAGE'],'.k',markersize=0.7)
    plt.title('NUV SExtractor, N={0}'.format(len(df)))
    plt.ylabel('FWHM [pix]')
    plt.ylim(0,8)
    plt.xlim(12,22)

    plt.suptitle(scan)
    plt.tight_layout()
    #plt.subplots_adjust(top=0.9)

    plt.savefig('/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/'+scan+'_match_new.png', dpi=190)
    plt.clf()

if __name__ == '__main__':
    if False:
        scans = ['0023', '0032', '0203', '0239', '0446', '0464', '0473', '0806', '0815', '1301', '1310', '1319',\
                '1616', '1679', '2192', '2714', '2750', '3236', '3245', '3281']
        ins = glob.glob('/scratch/dw1519/galex/fits/scan_map/11-27-2017/single/*_in.fits')
        date = '11-27-2017'
        ais = pyfits.open('../data/MyTable_4_jvc2688.fit')[1].data
        ais_cat = SkyCoord(ais['glon']*u.degree, ais['glat']*u.degree, frame='galactic')
        star_list = []
        for inm in ins:
            fname = os.path.basename(inm)
            scan = re.split('_', fname)[2]
            print(scan)
            cfile='/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
            df = load_data.load_catalog(cfile)
            c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
            print(c.shape)

            idx, d2d, d3d = c.match_to_catalog_sky(ais_cat)  
            mask=(d2d<0.001*u.degree) #& (df['gl']<50)
            N = np.sum(mask)

            gl = np.array(df['gl'][mask])-float(scan)/10.
            gb = np.array(df['gb'][mask])
            nuv = np.array(df['nuv'][mask])
            dif = df['nuv'][mask]-ais['nuv'][idx[mask]]
            cr_star = np.zeros(gl.shape[0])-1.

            print gl.shape
            print gb.shape
            print cr_star.shape
            count_mask = np.load('/scratch/dw1519/galex/data/star_photon/count_rate/'+scan+'_mask.npy')
            count_rate = np.load('/scratch/dw1519/galex/data/star_photon/count_rate/'+scan+'_ima.npy')

            for i in range(count_rate.shape[0]):
                cr = count_rate[i]
                print cr
                mask = count_mask[i]
                print mask
                gb_mask = (gb>=mask[0]) & (gb<=mask[1])
                print np.sum(gb_mask)
                #print gb_mask
                cr_star[gb_mask] = cr
            #print cr_star
            print np.sum(cr_star>0)
            star_list.append(np.column_stack([gl,gb,dif,nuv,cr_star]))
        star_list = np.concatenate(star_list, axis=0)
        np.save('/scratch/dw1519/galex/data/star_photon/count_rate/star_list.npy', star_list)

    if True:
        star_list = np.load('/scratch/dw1519/galex/data/star_photon/count_rate/star_list.npy')
        print star_list.shape
        mask = star_list[:,-1]>0.
        star_list = star_list[mask]
        print star_list.shape
        bins = np.arange(6)/100.
        bins[-1] = 1.
        bins = [0,0.015,0.02,0.023,0.026,0.03,0.05]
        plt.hist(star_list[:,-1], bins=bins)
        plt.show()
        
        flat0 = np.load('/scratch/dw1519/galex/data/star_photon/super_iterate3/flat.npy')
        date = '11-27-2017'
        median = np.median(star_list[:,2])
        fx = np.repeat(np.arange(800),800).reshape((800,800))
        fy = np.tile(np.arange(800),800).reshape((800,800))
        for i in range(6):
            flat = np.load('/scratch/dw1519/galex/data/star_photon/super_iterate{0}/flat9.npy'.format(i+10))
            r = np.sqrt((fx-399.5)**2+(fy-399.5)**2)
            flat[r>400] = 0
            data = np.load('/scratch/dw1519/galex/data/star_photon/super_iterate{0}/data.npy'.format(i+10))
            print 'data:{0}'.format(np.sum(data))
            profile0 = np.zeros(800)
            profile = np.zeros(800)
            n = 9.
            for rot in np.arange(n)*40.:
                flat0_tmp = rotate(flat0, rot,reshape=False,order=1,prefilter=False)
                flat_tmp = rotate(flat, rot,reshape=False,order=1,prefilter=False)
                profile0 += convert1(flat0_tmp)
                profile += convert1(flat_tmp)
            profile0 /= n
            profile /= n
            px = (np.arange(800)-399.5)*0.001666

            norm = np.sum(profile0[300:500])/np.sum(profile[300:500])
            print 'norm:{0}'.format(norm)
            flat = flat*norm
            np.save('/scratch/dw1519/galex/data/star_photon/super_iterate{0}/flat.npy'.format(i+10), flat)

            profile0 = np.zeros(800)
            profile = np.zeros(800)
            n = 9.
            for rot in np.arange(n)*40.:
                flat0_tmp = rotate(flat0, rot,reshape=False,order=1,prefilter=False)
                flat_tmp = rotate(flat, rot,reshape=False,order=1,prefilter=False)
                profile0 += convert1(flat0_tmp)
                profile += convert1(flat_tmp)
            profile0 /= n
            profile /= n

            cr_mask = (star_list[:,-1]>=bins[i]) & (star_list[:,-1]<bins[i+1])
            if np.sum(cr_mask) == 0 :
                continue

            f, ax1 = plt.subplots(1, sharex=True)
            ran=(-0.6,0.6)
            stars = star_list[cr_mask]
            dif = stars[:,2]
            gl = stars[:,0]
            #print gl.shape
            nuv = stars[:,3]
            #print nuv.shape
            #median= np.median(dif)
            #print median
            dif = dif-median
            dif_mask = (dif>=-0.5) & (dif<=0.5)
            dif = dif[dif_mask]#+median
            gl = gl[dif_mask]
            nuv = nuv[dif_mask]
            #print dif.shape
            #print gl.shape
            #print np.mean(dif), np.median(dif)
            ax1.scatter(gl, dif,\
                #c=nuv,\
                alpha=0.8, edgecolors='face', s=0.5, cmap=plt.get_cmap('jet'))
                # c=ais['glat'][idx[mask]], alpha=0.8, edgecolors='face', s=0.5, cmap=plt.get_cmap('jet'), vmin=-10, vmax=10)
            sm, b1 = np.histogram(gl, bins=120, range=ran, weights=dif)
            num, b2 = np.histogram(gl, bins=120, range=ran)

            x = np.convolve(b2, np.ones(2)/2., mode='valid')
            ma = sm/num
            ma = np.convolve(ma, np.ones(7)/7., mode='same')
            #print x
            #print ma
            ax1.plot(x, ma, '.r')
            ax1.axhline(y=0, color='k', linestyle='dashed', linewidth=1)

            ax1.plot(px, -2.5*np.log10(profile/profile0),'.k')

            ax1.set_ylim(-.75,.75)
            ax1.set_ylabel('NUV (SEx-AIS)')

            plt.xlim(-0.65,0.65)
            plt.title('average intensity:{0}-{1}'.format(bins[i], bins[i+1]))
            plt.xlabel('gl-scan_center')

            #plt.suptitle(scan)
            plt.tight_layout()
            f.subplots_adjust(hspace=0)

            plt.savefig('/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/cr_test_new{0}.png'.format(i), dpi=190)
            plt.clf()




