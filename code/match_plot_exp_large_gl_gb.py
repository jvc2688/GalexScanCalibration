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

        scans = ['0-10','9-19', '18-28', '27-37', '45-55' ]#['9-19']#['0050']#['0-10ds'] #['0050']
        #scans = ['72-82', '81-91', '90-100', '99-109','108-118', '117-127', '126-136', '135-145']
        #scans = ['144-154', '153-163', '162-172', '171-181','180-190']
        scans = ['198-208','207-217', '216-226', '225-235', '234-244' ]
        scans = ['252-262','261-271', '270-280', '279-289', '288-298', '297-307' ]
        #scans = ['315-325', '324-334', '333-343', '351-1' ]
        scans = ['0-10', '18-28', '45-55', '72-82']#, '90-100', '135-145', '198-208', '216-226', '234-244', '297-307']
        scans = ['261-271']
        #scans = ['1634', '2183', '2552', '2561', '2570', '2615', '2642', '2921', '3497']
        #scans = ['2174']
        ranges = [(0.,10.), (9.,19.), (18.,28.), (27.,37.), (45,55)]
        #ranges = [(72.,82.), (81.,91.), (90.,100.), (99.,100.), (108,118), (117,127), (126,136), (135,145)]
        #ranges = [(144.,154.), (153.,163.), (162.,172.), (171.,181.), (180,190)]
        ranges = [(198.,208.), (207.,217.), (216.,226.)]
        ranges = [(252.,262.), (261.,271.), (270.,280.), (279.,289.), (288,298), (297,307), (225.,235.), (234,244)]
        #ranges = [(315.,325.), (324.,334.), (333.,343.), (351.,1.)]


        ran = (234,244)
        date = '12-06-2017'
        nrange = 0
        nuv_sex_list =[]
        nuv_ais_list = []
        gl_list = []
        gb_list = []
        scan = '261-271'

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

        gl = df['gl']
        gb = df['gb']
        #dif = df['nuv'][mask]-ais['nuv'][idx[mask]]

        f, (ax1) = plt.subplots(1, sharex=True, figsize=(9.5,19.5))
        cax = ax1.scatter(gl, gb,
            alpha=1., edgecolors='face', s=0.5, cmap=plt.get_cmap('jet'))


        plt.tight_layout()
        f.subplots_adjust(hspace=0)

        plt.savefig('/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/261-271_gl_gb.png', dpi=190)
        plt.clf()



