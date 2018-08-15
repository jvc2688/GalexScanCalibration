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
      q_mask = data[:,3]>5
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
    return np.array([data[:,0],ra,dec, data[:,4], data[:,5]]).T

def detector2gon(pos_list, size):
  pos_array = np.array(pos_list, dtype='float64')
  pos = (pos_array/50.-0.5)*36000*800*0.001666
  return pos


if __name__ == '__main__':

    scans = ['0050']
    date = '08-21-2017'#'08-21-2017'
    for scan in scans:
        print(scan)
        cfile='/scratch/dw1519/galex/fits/scan_map/catalog/'+date+'/starcat_'+scan+'mapweight_fec_fwhm.txt'
        df = load_data.load_catalog(cfile)
        #nuv_mask = (df['nuv'] >= 10.0) & (df['nuv']<= 17.0)
        #df = df[nuv_mask]
        c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
        c_radec = c.transform_to('fk5')
        print(c.shape)

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

            csvs = glob.glob("{0}/{1}/split/*.csv".format(data_path, name))
            photon_list = []
            for f in csvs:
                print(f)
                data = np.loadtxt(f,delimiter=',', usecols=[0,3,4,5,6,7])
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
                idxc, idxcatalog, d2d_p, d3d_p = sky_data.search_around_sky(c, 0.0014*u.deg)
                if len(idxc)>0:
                    photon_list.append(np.column_stack([idxc, data[idxcatalog,0], data[idxcatalog,3:]]))
                else:
                    print 'no star'
            photon_list = np.concatenate(photon_list, axis=0)
        np.save('/scratch/dw1519/galex/data/star_photon/'+scan+'_photon.npy', photon_list)
        np.save('/scratch/dw1519/galex/data/star_photon/'+scan+'_star.npy', np.array([df['gl'], df['gb'], df['nuv']]).T)


