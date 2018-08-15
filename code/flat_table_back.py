import numpy as np
import matplotlib.pyplot as plt
import sys
from fast_histogram import histogram2d
from scipy import interpolate
import gnomonic as gn
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
import glob
import load_data
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits as pyfits


def dis_correct(data, asp, dis_map=None, ya_corr=None, q_corr=None, cut=False):
    if cut:
      ya_mask = data[:,2]>=2
      q_mask = data[:,3]>=0
      data = data[ya_mask&q_mask]
    ix = np.digitize(data[:,0]/1000., asp[:,0])-1

    if dis_map is not None:
      low_mask = data[:,1]<16
      high_mask = data[:,1]>=16

      coo = data[low_mask, 4:6]
      data[low_mask, 4] -= dis_map['xi_l'](coo)
      data[low_mask, 5] -= dis_map['eta_l'](coo)

      coo = data[high_mask, 4:6]
      data[high_mask, 4] -= dis_map['xi_h'](coo)
      data[high_mask, 5] -= dis_map['eta_h'](coo)

    if q_corr is not None:
      q_mask = data[:,3]<24
      data[q_mask,4] -= q_corr[data[q_mask,3].astype(int),-2]*36000*800*0.001666/2400.
      data[q_mask,5] -= q_corr[data[q_mask,3].astype(int),-1]*36000*800*0.001666/2400.

    if ya_corr is not None:
      data[:,4] -= ya_corr[data[:,2].astype(int)-2,-2]*36000*800*0.001666/2400.
      data[:,5] -= ya_corr[data[:,2].astype(int)-2,-1]*36000*800*0.001666/2400.

    ra, dec = gn.gnomrev_simple(data[:,4], data[:,5], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return np.array([data[:,0],ra,dec, data[:,4], data[:,5], data[:,3]]).T

def detector2gon(pos_list, size):
  pos_array = np.array(pos_list, dtype='float64')
  pos = (pos_array/50.-0.5)*36000*800*0.001666
  return pos


if __name__ == '__main__':

    scans = ['2120', '2633','3047']
    scans = ['0158', '0689','1247','1508', '1778', '3326', '3587']
    scans = ['0212','0302','0500','0743','2291','2453', '2759', '2885', '2930', '3155', '3227', '3497']
    scans = ['0005', '0023', '0446', '0464','0482', '0806', '0941', '0995']
    scans = ['0023', '0212', '0302', '0464', '0500', '0689', '0806', '0941', '1247', '1508', '1778', 
            '2120', '2291', '2453','2759', '2930', '3155', '3227', '3326', '3497']
    scans = ['0023', '0212']
    scans = ['0023', '0212', '0464','0689']
    date = '08-21-2017'#'08-21-2017'
    scan_f = sys.argv[1]
    scans = []
    with open(scan_f) as f:
      scans = f.read().splitlines()
    print scans
    size=800
    Nb = 6000
    ap = 0.004

    for scan in scans:
        print(scan)
        cfile='/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
        df = load_data.load_catalog(cfile)
        print len(df)
        mask = df['nuv']<14.5
        c1 = SkyCoord(df['gl'][mask]*u.degree, df['gb'][mask]*u.degree, frame='galactic')
        mask = (df['nuv']>=14.5) & (df['nuv']<16.5)
        c2 = SkyCoord(df['gl'][mask]*u.degree, df['gb'][mask]*u.degree, frame='galactic')
        mask = (df['nuv']>=16.5)
        c3 = SkyCoord(df['gl'][mask]*u.degree, df['gb'][mask]*u.degree, frame='galactic')        
        print(c1.shape[0]+c2.shape[0]+c3.shape[0])

        gl = np.array(df['gl'])
        gb = np.array(df['gb'])
        bmax = np.max(gb)
        bmin = np.min(gb)
        lmax = np.max(gl)
        lmin = np.min(gl)
        if lmax-lmin > 2:
            gl[gl<350] = 360.-gl[gl<350]
            lmax = np.max(gl)
            lmin = np.min(gl)

        bl = np.random.uniform(lmin+0.02, lmax-0.02, Nb)
        bb = np.random.uniform(bmin+1.,bmax-1.,Nb)
        bc = SkyCoord(bl*u.degree, bb*u.degree, frame='galactic')

        bmask = np.ones(bc.shape[0], dtype=bool)
        idxc1, idxcatalog1, d2d_p1, d3d_p1 = c1.search_around_sky(bc, (0.02+ap)*u.deg)
        idxc2, idxcatalog2, d2d_p2, d3d_p2 = c2.search_around_sky(bc, (0.015+ap)*u.deg)
        idxc3, idxcatalog3, d2d_p3, d3d_p3 = c3.search_around_sky(bc, (0.01+ap)*u.deg)

        bmask[idxc1] = False
        bmask[idxc2] = False
        bmask[idxc3] = False

        bc = bc[bmask]
        np.save('/scratch/dw1519/galex/data/star_photon/back/'+scan+'_bstar.npy', np.array([bc.l.deg, bc.b.deg]).T)
        print bc.shape

        with open('../name_scan/%s'%scan) as f:
          name_list = f.read().splitlines()
        print name_list

        data_path = '/beegfs/dw1519/galex/data'
        asprta = '/scratch/dw1519/galex/AIS_GAL_SCAN/asprta'
        asprta_suffix = '-cal-sec-dis-cut-correct-nuv-new'
        resolution = 0.5
        scst = '/scratch/dw1519/galex/AIS_GAL_SCAN/scst'
        ya = '/scratch/dw1519/galex/data'
        q = '/scratch/dw1519/galex/data'
        cut = True

        num= int(scan)
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
        #    content = re.sub(r'^dis=.*$', 'dis=/scratch/dw1519/galex/data/distortion/896-1157', content, flags=re.MULTILINE)
        elif num>=1166 and num<=1373:
            dis_l = '/scratch/dw1519/galex/data/distortion/1166-1373-xa_low-centroids.npy'
            dis_h = '/scratch/dw1519/galex/data/distortion/1166-1373-xa_high-centroids.npy'        
        else:
            dis_l = '/scratch/dw1519/galex/data/distortion/194-437-xa_low-centroids.npy'
            dis_h = '/scratch/dw1519/galex/data/distortion/194-437-xa_high-centroids.npy' 

        if dis_l is not None:
            detector_size = [50,50]
            dis_map_l = np.load(dis_l)
            dis_map_h = np.load(dis_h)
            dis_map_l = dis_map_l.reshape(50*50,2)/detector_size[0]*36000*800*0.001666
            dis_map_h = dis_map_h.reshape(50*50,2)/detector_size[0]*36000*800*0.001666
            cen_x = np.arange(detector_size[0])+0.5
            cen_y = np.arange(detector_size[1])+0.5
            xi = np.repeat(detector2gon(cen_x, detector_size[0]), detector_size[1])
            eta = np.tile(detector2gon(cen_y, detector_size[1]), detector_size[0])
            points = np.array([xi,eta]).T

            fxi_l = interpolate.LinearNDInterpolator(points, dis_map_l[:,0])
            feta_l = interpolate.LinearNDInterpolator(points, dis_map_l[:,1])
            fxi_h = interpolate.LinearNDInterpolator(points, dis_map_h[:,0])
            feta_h = interpolate.LinearNDInterpolator(points, dis_map_h[:,1])

            dis_map = {'xi_l':fxi_l, 'xi_h':fxi_h, 'eta_l':feta_l, 'eta_h':feta_h}
        else:
            dis_map = None

        for name in name_list:
            co_data = pyfits.open('{0}/{1}{2}-asprta.fits'.format(asprta, name, asprta_suffix))[1].data
            n = resolution/0.005
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

            scst_file = pyfits.open('../AIS_GAL_SCAN/scst/{0}-scst.fits'.format(name))
            scst_data = scst_file[1].data
            scst_time = scst_data['pktime'].copy()
            hv = scst_data['hvnom_nuv'].copy()
            scst_file.close()

            if q is not None:
                q_corr = np.load(q+'/'+name+'/q_corr.npy')
            else:
                q_corr = None
            if ya is not None:
                ya_corr = np.load(ya+'/'+name+'/ya_corr.npy')
            else:
                ya_corr = None

            csvs = glob.glob("{0}/{1}/split/*.csv".format(data_path, name))
            photon_list = []
            for f in csvs:
                print(f)
                data = np.loadtxt(f, delimiter=',', usecols=[0,3,4,5,6,7])
                if len(data.shape)<2:
                    continue
                scst_ix = np.digitize(data[:,0]/1000., scst_time)
                ix_mask = scst_ix<scst_time.shape[0]
                scst_ix = scst_ix[ix_mask]
                data = data[ix_mask]
                data = data[hv[scst_ix]>0]
                if len(data)==0:
                    print('skip')
                    continue
                data = dis_correct(data, asp, dis_map, ya_corr, q_corr, cut)
                sky_data = SkyCoord(data[:,1:3], unit='deg', frame='fk5', equinox='J2000.0')
                #gal = sky_data.transform_to('galactic')
                idxc, idxcatalog, d2d_p, d3d_p = sky_data.search_around_sky(bc, ap*u.deg)
                if len(idxc)>0:
                    photon_list.append(np.column_stack([idxc, data[idxcatalog,0], data[idxcatalog,3:]]))
                else:
                    print 'no star'
            photon_list = np.concatenate(photon_list, axis=0)
            np.save('/scratch/dw1519/galex/data/star_photon/back/'+name+'_photon.npy', photon_list)
            #np.save('/scratch/dw1519/galex/data/star_photon/back/'+name+'_bstar.npy', np.array([df['gl'], df['gb'], df['nuv']]).T)


