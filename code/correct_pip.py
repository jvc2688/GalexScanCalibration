import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import c3
from scipy import stats
import glob
import re
import os
import load_data as ld
from scipy.interpolate import splev, splrep
import gnomonic as gn
import gc
from astropy.io import fits as pyfits
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys

def get_key(path):
    return int(re.split('_|/|\.', path)[-2])

def get_key_info(path):
    return int(re.split('_|\.', path)[-2])


def correct(name):
    #name = sys.argv[1]
    print name
    scan = '{0:04d}'.format(int(re.split('_', name)[3])) #'3191'
    print scan
    #name = 'AIS_GAL_SCAN_0'+scan+'_0002'
    out = '/scratch/dw1519/galex/data/'+name
    date = '08-21-2017_1'
    cfile = '/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
    asprta = '/scratch/dw1519/galex/AIS_GAL_SCAN/asprta/'+name+'-cal-sec-dis-cut-asprta.fits'
    scst = '/scratch/dw1519/galex/AIS_GAL_SCAN/scst/'+name+'-scst.fits'
    num = int(scan)
    if num>=5 and num<=185:
        dis_l = '/scratch/dw1519/galex/data/distortion/5-185-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/5-185-xa_high-centroids.npy'
    elif num>=194 and num<=437:
        dis_l = '/scratch/dw1519/galex/data/distortion/194-437-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/194-437-xa_high-centroids.npy'
    elif num>=446 and num<=887:
        dis_l = '/scratch/dw1519/galex/data/distortion/446-887-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/446-887-xa_high-centroids.npy'
    #elif num>=896 and num<=1157:
    #    dis_l = '/scratch/dw1519/galex/data/distortion/896-1157-xa_low-centroids.npy'
    #    dis_h = '/scratch/dw1519/galex/data/distortion/896-1157-xa_high-centroids.npy'
    elif num>=1166 and num<=1373:
        dis_l = '/scratch/dw1519/galex/data/distortion/1166-1373-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/1166-1373-xa_high-centroids.npy'
    else:
        dis_l = '/scratch/dw1519/galex/data/distortion/194-437-xa_low-centroids.npy'
        dis_h = '/scratch/dw1519/galex/data/distortion/194-437-xa_high-centroids.npy'
    print dis_l
    print dis_h
    data = np.load('/scratch/dw1519/galex/data/photons/'+name+'.npy')
    print('data:{0}'.format(data.shape))
    tycho2 = pyfits.open('/scratch/dw1519/galex/data/tycho2.fits')[1].data
    df = ld.load_catalog(cfile)
    c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
    catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask=d2d<0.001*u.degree

    '''
    dtype = np.dtype([('tycho_num', int), ('Glon', '>f4'), ('Glat', '>f4'), ('RAJ2000', '>f4'), ('DEJ2000', '>f4'),
             ('flux', float), ('nuv', float), ('gl', float), ('gb', float)])
    matched_tycho = tycho2[idx[mask]]
    matched_df = df[mask]
    matched_catalog = np.core.records.fromarrays(np.array([idx[mask], matched_tycho['Glon'], matched_tycho['Glat'],
            matched_tycho['RAJ2000'], matched_tycho['DEJ2000'], np.array(matched_df['FLUX_AUTO']),
            np.array(matched_df['nuv']), np.array(matched_df['gl']),
            np.array(matched_df['gb'])]), dtype=dtype)
    '''

    resolution = 0.5
    hdulist = pyfits.open(asprta)
    co_data = hdulist[1].data
    n=resolution/0.005 #100
    T = co_data['T']
    ra = co_data['ra']
    dec = co_data['dec']
    roll = co_data['roll']
    t_new = np.arange((T.shape[0]-1)*n)*0.005+T[0]
    spl_ra = splrep(T, ra)
    spl_dec = splrep(T, dec)
    spl_roll = splrep(T, roll)
    ra_new = splev(t_new, spl_ra)
    dec_new = splev(t_new, spl_dec)
    roll_new = splev(t_new, spl_roll)

    asp = np.array([t_new, ra_new, dec_new, roll_new]).T
    #np.save('../data/photon_list/%s_asp_fix.npy'%name, asp_cal)
    hdulist.close()
    
    #asp = np.load('../data/photon_list/{0}{1}-correct_asp_new.npy'.format(name, suffix))
    scst_file = pyfits.open(scst)
    scst_data = scst_file[1].data
    scst_time = scst_data['pktime'].copy()
    hv = scst_data['hvnom_nuv'].copy()
    scst_file.close()

    dis_map = ld.load_dis(dis_l, dis_h)

    scst_ix = np.digitize(data[:,0]/1000., scst_time)
    ix_mask = scst_ix<scst_time.shape[0]
    scst_ix = scst_ix[ix_mask]
    data = data[ix_mask]
    data = data[hv[scst_ix]>0]
    
    dx_list = []
    dy_list = []
    ya_list = []
    q_list = []
    n = 50
    step = int(data.shape[0]/n)
    print step
    for i in range(n-1):
        print i
        data_tmp = ld.dis_correct(data[i*step:(i+1)*step], dis_map, asp, None, None)
        sky_data = SkyCoord(data_tmp[:,0:2], unit='deg', frame='fk5', equinox='J2000.0')
        print 'sky_data: {0}'.format(sky_data.shape)
        idxc, idxcatalog, d2d_p, d3d_p = sky_data.search_around_sky(catalog[idx[mask]], 0.0045*u.deg)
        print 'idxc: {0}'.format(idxc.shape)
        xi, eta = gn.gnomfwd_simple(tycho2['RAJ2000'][idx[mask][idxc]], tycho2['DEJ2000'][idx[mask][idxc]],
                 data_tmp[idxcatalog,-3], data_tmp[idxcatalog,-2], -data_tmp[idxcatalog,-1], 1/36000., 0)

        dx = (data_tmp[idxcatalog,6]-xi)/36000./800/0.001666*2400
        dy = (data_tmp[idxcatalog,7]-eta)/36000./800/0.001666*2400
        ya = data_tmp[idxcatalog,4]
        q = data_tmp[idxcatalog,5]

        dx_list.append(dx)
        dy_list.append(dy)
        ya_list.append(ya)
        q_list.append(q)
        data_tmp=sky_data=xi=eta=None
        gc.collect()
    
    print n-1
    data_tmp = ld.dis_correct(data[(n-1)*step:], dis_map, asp, None, None)
    sky_data = SkyCoord(data_tmp[:,0:2], unit='deg', frame='fk5', equinox='J2000.0')
    print 'sky_data: {0}'.format(sky_data.shape)
    idxc, idxcatalog, d2d_p, d3d_p = sky_data.search_around_sky(catalog[idx[mask]], 0.0045*u.deg)
    print 'idxc: {0}'.format(idxc.shape)
    xi, eta = gn.gnomfwd_simple(tycho2['RAJ2000'][idx[mask][idxc]], tycho2['DEJ2000'][idx[mask][idxc]],
             data_tmp[idxcatalog,-3], data_tmp[idxcatalog,-2], -data_tmp[idxcatalog,-1], 1/36000., 0)

    dx = (data_tmp[idxcatalog,6]-xi)/36000./800/0.001666*2400
    dy = (data_tmp[idxcatalog,7]-eta)/36000./800/0.001666*2400
    ya = data_tmp[idxcatalog,4]
    q = data_tmp[idxcatalog,5]

    dx_list.append(dx)
    dy_list.append(dy)
    ya_list.append(ya)
    q_list.append(q)
    data_tmp=sky_data=xi=eta=None
    gc.collect()

    dx = np.concatenate(dx_list, axis=0)
    dy = np.concatenate(dy_list, axis=0)
    ya = np.concatenate(ya_list, axis=0)
    q = np.concatenate(q_list, axis=0)

    data=scst_ix=ix_mask=None
    gc.collect()
    print dx.shape
    print 'max dy: {0}'.format(np.max(dy))
    print 'max dx: {0}'.format(np.max(dx))
    mask = (dy>=-8) & (dy<=8) & (q<24)
    imsz = [24, 64]
    H,xedges,yedges=np.histogram2d(q[mask], dy[mask],\
                             bins=imsz)

    plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
                 extent=[-8, 8, 0, 24], origin='lower', norm=LogNorm())

    dy_q = []
    for i in range(24):
        q_mask = (q[mask] == i)
        dy_c = np.median(dy[mask][q_mask])
        dy_q.append(dy_c)
    plt.plot(dy_q, np.arange(24), '.r')
    print('dy_q:')
    print dy_q

    plt.xlim(-8,8)
    plt.ylim(0,24)
    plt.title('scan{0:>5}'.format(scan)) 
    plt.xlabel(r'$\Delta y$ [pixel]')
    plt.ylabel('q')
    plt.savefig(out+'/q-y_tot.pdf', dpi=190)
    plt.clf()

    plt.plot(np.arange(24), dy_q, '.r')
    plt.savefig(out+'/q-y_c.pdf', dpi=190)
    plt.clf()

    mask = (dx>=-8) & (dx<=8) & (q<24)
    imsz = [24, 64]
    H,xedges,yedges=np.histogram2d(q[mask], dx[mask],\
                             bins=imsz)

    plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
                 extent=[-8, 8, 0, 24], origin='lower', norm=LogNorm())


    dx_q = []
    for i in range(24):
        q_mask = (q[mask] == i)
        dx_c = np.median(dx[mask][q_mask])
        dx_q.append(dx_c)
    print('dx_q:')
    print dx_q
    plt.plot(dx_q, np.arange(24), '.r')

    plt.xlim(-8,8)
    plt.ylim(0,24)
    plt.title('scan{0:>5}'.format(scan)) 
    plt.xlabel(r'$\Delta x$ [pixel]')
    plt.ylabel('q')
    plt.savefig(out+'/q-x_tot.pdf', dpi=190)
    plt.clf()

    plt.plot(np.arange(24), dx_q, '.r')
    plt.savefig(out+'/q-x_c.pdf', dpi=190)
    plt.clf()

    q_corr = np.array([np.arange(24), np.array(dx_q), np.array(dy_q)]).T
    np.save(out+'/q_corr.npy', q_corr)
    q_mask = q<24
    dy[q_mask] = dy[q_mask]-q_corr[q[q_mask].astype(int),-1]
    dx[q_mask] = dx[q_mask]-q_corr[q[q_mask].astype(int),-2]

    mask = (dy>=-8) & (dy<=8)
    imsz = [30, 64]
    H,xedges,yedges=np.histogram2d(ya[mask], dy[mask],\
                             bins=imsz)

    plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
                 extent=[-8, 8, 2, 32], origin='lower', norm=LogNorm())

    dy_ya = []
    for i in range(2,32):
        ya_mask = (ya[mask] == i)
        dy_c = np.median(dy[mask][ya_mask])
        dy_ya.append(dy_c)
    print('dy_ya:')
    print dy_ya
    plt.plot(dy_ya, np.arange(30)+2, '.r')

    plt.xlim(-8,8)
    plt.ylim(2,32)
    plt.title('scan{0:>5}'.format(scan)) 
    plt.xlabel(r'$\Delta y$ [pixel]')
    plt.ylabel('ya')
    plt.savefig(out+'/ya-y_tot.pdf', dpi=190)
    plt.clf()

    plt.plot(np.arange(30)+2, dy_ya, '.k')
    plt.savefig(out+'/ya-y_c.pdf', dpi=190)
    plt.clf()

    mask = (dx>=-8) & (dx<=8)
    imsz = [30, 64]
    H,xedges,yedges=np.histogram2d(ya[mask], dx[mask],\
                             bins=imsz)

    plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
                 extent=[-8, 8, 2, 32], origin='lower', norm=LogNorm())

    dx_ya = []
    for i in range(2,32):
        ya_mask = (ya[mask] == i)
        dx_c = np.median(dx[mask][ya_mask])
        dx_ya.append(dx_c)
    print('dx_ya:')
    print dx_ya
    plt.plot(dx_ya, np.arange(30)+2, '.r')

    plt.xlim(-8,8)
    plt.ylim(2,32)
    plt.title('scan{0:>5}'.format(scan)) 
    plt.xlabel(r'$\Delta x$ [pixel]')
    plt.ylabel('ya')
    plt.savefig(out+'/ya-x_tot.pdf', dpi=190)
    plt.clf()

    plt.plot(np.arange(30)+2, np.array(dx_ya), '.k')
    plt.savefig(out+'/ya-x_c.pdf', dpi=190)
    plt.clf()

    ya_corr = np.array([np.arange(30)+2, np.array(dx_ya), np.array(dy_ya)]).T
    np.save(out+'/ya_corr.npy', ya_corr)
    ya_mask = q<24
    dy = dy-ya_corr[ya.astype(int)-2,-1]
    dx = dx-ya_corr[ya.astype(int)-2,-2]

    #plot corrected
    mask = (dy>=-8) & (dy<=8) & (q<24)
    imsz = [24, 64]
    H,xedges,yedges=np.histogram2d(q[mask], dy[mask],\
                             bins=imsz)

    plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
                 extent=[-8, 8, 0, 24], origin='lower', norm=LogNorm())

    dy_cs = []
    for i in range(24):
        q_mask = q[mask] == i
        dy_c = np.median(dy[mask][q_mask])
        dy_cs.append(dy_c)
    print('dy_q:')
    print dy_cs
    plt.plot(dy_cs, np.arange(24), '.r')

    plt.xlim(-8,8)
    plt.ylim(0,24)
    plt.title('scan{0:>5}'.format(scan)) 
    plt.xlabel(r'$\Delta y$ [pixel]')
    plt.ylabel('q')
    plt.savefig(out+'/q-y_tot_new.pdf', dpi=190)
    plt.clf()
    #np.save(out+'/q-y.npy', np.array([np.arange(24),dy_cs]).T)

    plt.plot(np.arange(24), dy_cs, '.r')
    plt.savefig(out+'/q-y_c_new.pdf', dpi=190)
    plt.clf()

    mask = (dx>=-8) & (dx<=8) & (q<24)
    imsz = [24, 64]
    H,xedges,yedges=np.histogram2d(q[mask], dx[mask],\
                             bins=imsz)

    plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
                 extent=[-8, 8, 0, 24], origin='lower', norm=LogNorm())


    dy_cs = []
    for i in range(24):
        q_mask = q[mask] == i
        dy_c = np.median(dx[mask][q_mask])
        dy_cs.append(dy_c)
    print('dx_q:')
    print dy_cs
    plt.plot(dy_cs, np.arange(24), '.r')

    plt.xlim(-8,8)
    plt.ylim(0,24)
    plt.title('scan{0:>5}'.format(scan)) 
    plt.xlabel(r'$\Delta x$ [pixel]')
    plt.ylabel('q')
    plt.savefig(out+'/q-x_tot_new.pdf', dpi=190)
    plt.clf()
    #np.save(out+'/q-x.npy', np.array([np.arange(24),dy_cs]).T)

    plt.plot(np.arange(24), dy_cs, '.r')
    plt.savefig(out+'/q-x_c_new.pdf', dpi=190)
    plt.clf()


    mask = (dy>=-8) & (dy<=8)
    imsz = [30, 64]
    H,xedges,yedges=np.histogram2d(ya[mask], dy[mask],\
                             bins=imsz)

    plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
                 extent=[-8, 8, 2, 32], origin='lower', norm=LogNorm())

    dy_cs = []
    for i in range(2,32):
        ya_mask = ya[mask] == i
        dy_c = np.median(dy[mask][ya_mask])
        dy_cs.append(dy_c)
    print('dy_ya:')
    print dy_cs
    plt.plot(dy_cs, np.arange(30)+2, '.r')

    plt.xlim(-8,8)
    plt.ylim(2,32)
    plt.title('scan{0:>5}'.format(scan)) 
    plt.xlabel(r'$\Delta y$ [pixel]')
    plt.ylabel('ya')
    plt.savefig(out+'/ya-y_tot_new.pdf', dpi=190)
    plt.clf()

    plt.plot(np.arange(30)+2, dy_cs, '.k')
    plt.savefig(out+'/ya-y_c_new.pdf', dpi=190)
    plt.clf()
    #np.save(out+'/ya-y.npy', np.array([np.arange(30)+2,dy_cs]).T)

    mask = (dx>=-8) & (dx<=8)
    imsz = [30, 64]
    H,xedges,yedges=np.histogram2d(ya[mask], dx[mask],\
                             bins=imsz)

    plt.imshow(H,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal',\
                 extent=[-8, 8, 2, 32], origin='lower', norm=LogNorm())

    dy_cs = []
    for i in range(2,32):
        ya_mask = ya[mask] == i
        dy_c = np.median(dx[mask][ya_mask])
        dy_cs.append(dy_c)
    print('dx_ya:')
    print dy_cs
    plt.plot(dy_cs, np.arange(30)+2, '.r')

    plt.xlim(-8,8)
    plt.ylim(2,32)
    plt.title('scan{0:>5}'.format(scan)) 
    plt.xlabel(r'$\Delta x$ [pixel]')
    plt.ylabel('ya')
    plt.savefig(out+'/ya-x_tot_new.pdf', dpi=190)
    plt.clf()

    plt.plot(np.arange(30)+2, dy_cs, '.k')
    plt.savefig(out+'/ya-x_c_new.pdf', dpi=190)
    plt.clf()

if __name__ == '__main__':


    if True:
        sub_name = sys.argv[1]
        print sub_name
        date = '08-21-2017_1'
        scan_list = []
        with open(sub_name, 'r') as f:
            scan_list = f.read().splitlines()
        for name in scan_list:
            print name
            num = int(re.split('_', name)[3])
            scan = '%04d'%num
            print scan
            catalog = '/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
            ya = '/scratch/dw1519/galex/data/'+name+'/ya_corr.npy'
            if os.path.isfile(catalog):
                if(not os.path.isfile(ya)):
                    correct(name)
                else:
                    print "existed"
            else:
                print 'no source file'
