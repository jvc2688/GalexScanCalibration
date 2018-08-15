from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import load_data
from matplotlib.colors import LogNorm
from astropy import wcs
#plt.style.use('seaborn-paper')

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
    if True:
        scans = ['0023', '0032', '0203', '0239', '0446', '0464', '0473', '0806', '0815', '1301', '1310', '1319',\
                '1616', '1679', '2192', '2714', '2750', '3236', '3245', '3281']

        scans = ['0050']#['0050']#['0-10ds'] #['0050']
        #scans = ['1634', '2183', '2552', '2561', '2570', '2615', '2642', '2921', '3497']
        #scans = ['2174']
        ran = (4.3,5.7)
        date = '11-21-2017'
        for scan in scans:
            print(scan)
            tycho2 = pyfits.open('../data/tycho2.fits')[1].data
            cfile='/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
            df = load_data.load_catalog(cfile)
            c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
            print(c.shape)
            catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
            idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
            mask=d2d<0.001*u.degree
            N = np.sum(mask)
            print(N)

            exp = pyfits.open('/scratch/dw1519/galex/fits/scan_map/'+date+'/count_map_'+scan+'_exp.fits')[0]
            w = wcs.WCS(exp.header)
            e = exp.data
            mask = e>0
            nmask = np.sum(mask, axis=0)
            nsort = np.argsort(nmask)
            sortm = nmask[nsort]
            n = nsort[sortm>10000][0]
            print nmask[n]
            e_mean = np.mean(e[e[:,n]>0], axis=0)

            x = np.arange(e.shape[1])
            y = np.zeros(e.shape[1])+e.shape[0]/2.
            egl, egb = w.all_pix2world(x,y,0)


            ais = pyfits.open('../data/MyTable_4_jvc2688.fit')[1].data
            ais_cat = SkyCoord(ais['glon']*u.degree, ais['glat']*u.degree, frame='galactic')
            idx, d2d, d3d = c.match_to_catalog_sky(ais_cat)  
            mask=(d2d<0.001*u.degree) #& (df['gl']<50)
            N = np.sum(mask)
            print('ais match:{0}'.format(N))
            #plt.subplot(2, 2, 2)
            gl = df['gl'][mask]
            gb = df['gb'][mask]
            dif = df['nuv'][mask]-ais['nuv'][idx[mask]]

            group = []
            g1 = (gb>-10) & (gb<=-5)
            group.append(g1)
            g2 = (gb>-5) & (gb<=0)
            group.append(g2)
            g3 = (gb>0) & (gb<=5)
            group.append(g3)
            g4 = (gb>5) & (gb<=10)
            group.append(g4)

            i = 0
            for g in group:
                f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(8,6))
                glg = gl[g]
                difg = dif[g]
                difg = difg-np.median(difg)
                dif_mask = (difg>=-0.5) & (difg<=0.5)
                difg = difg[dif_mask]
                glg = glg[dif_mask]
                print difg.shape
                print glg.shape
                print np.mean(difg), np.median(difg)
                ax1.scatter(glg,\
                    difg,\
                    #c=ais['nuv'][idx[mask]],\
                    alpha=0.8, edgecolors='face', s=0.5, cmap=plt.get_cmap('jet'), vmin=13, vmax=19)
                    # c=ais['glat'][idx[mask]], alpha=0.8, edgecolors='face', s=0.5, cmap=plt.get_cmap('jet'), vmin=-10, vmax=10)
                sm, b1 = np.histogram(glg, bins=500, range=ran, weights=difg)
                num, b2 = np.histogram(glg, bins=500, range=ran)

                x = np.convolve(b2, np.ones(2)/2., mode='valid')
                ma = sm/num
                ma = np.convolve(ma, np.ones(7)/7., mode='same')
                #print x
                #print ma
                ax1.plot(x, ma, '.r')
                ax1.axhline(y=0, color='k', linestyle='dashed', linewidth=1)

                ax1.set_ylim(-1.,1.)
                ax1.set_ylabel('NUV (SEx-AIS)')
                #plt.title('NUV AIS, N={0}'.format(N))

                ax2.scatter(egl,\
                    e_mean,\
                    #c=ais['nuv'][idx[mask]],\
                    alpha=0.8, edgecolors='face', s=0.5, cmap=plt.get_cmap('jet'), vmin=13, vmax=19)
                ax2.set_ylabel('EXP [s]')

                #nuv_list = [[14,15],[15,16],[16,17],[17,18],[18,19],[19,20]]
                #mask_list = [(e[0]<=ais['nuv'][idx[mask]]) & (e[1]>ais['nuv'][idx[mask]]) for e in nuv_list]
                #x = [np.median(ais['nuv'][idx[mask]][m]) for m in mask_list]
                #y = [np.median((df['nuv'][mask]-ais['nuv'][idx[mask]])[m]) for m in mask_list]
                #ye = [np.std((df['nuv'][mask]-ais['nuv'][idx[mask]])[m]) for m in mask_list]
                #plt.errorbar(x, y, yerr=ye)
                #angle = np.load('/scratch/dw1519/galex/fits/scan_map/'+date+'/'+scan+'.npy')
                #ax3.set_ylabel('Roll [deg]')
                #ax3.plot(angle[:,0], angle[:,1],'.')

                plt.xlim(4.3,5.7)
                #plt.xlim(28,35)
                #plt.xlim(172,181)
                #plt.xlim(0,7)
                #plt.ylim(-2,2)
                #plt.colorbar()
                plt.xlabel('gl')
                plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

                #plt.suptitle(scan)
                plt.tight_layout()
                f.subplots_adjust(hspace=0)
                #plt.subplots_adjust(top=0.9)

                plt.savefig('/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/'+scan+'_match_exp_new{0}.png'.format(i), dpi=190)
                plt.clf()
                i += 1
                '''
                np.save('/scratch/dw1519/galex/data/star_photon/cut_iterate6/stars_in.npy', 
                    np.column_stack([gl,df['nuv'][mask]-ais['nuv'][idx[mask]], ais['nuv'][idx[mask]], df['FWHM_IMAGE'][mask]]))
                np.save('/scratch/dw1519/galex/data/star_photon/cut_iterate6/stars_in_all.npy',
                    np.column_stack([df['gl'], df['nuv'], df['FWHM_IMAGE']]))
                '''



