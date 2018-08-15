import numpy as np
import matplotlib.pyplot as plt
import sys
#from fast_histogram import histogram2d
from scipy import interpolate
import gnomonic as gn
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
import glob
import load_data
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits as pyfits
from astropy.nddata.utils import block_reduce


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
    size=800

    f_scans = ['AIS_GAL_SCAN_00005_0001', 'AIS_GAL_SCAN_00446_0001', 'AIS_GAL_SCAN_00995_0001',
                'AIS_GAL_SCAN_00005_0001', 'AIS_GAL_SCAN_01247_0001', 'AIS_GAL_SCAN_02120_0002',
                'AIS_GAL_SCAN_03326_0002']
    full = np.zeros((size,size))
    cut = np.zeros((size,size))
    for scan in f_scans:
        full += np.load('/beegfs/dw1519/galex/data/'+scan+'/full.npy')
        cut += np.load('/beegfs/dw1519/galex/data/'+scan+'/cut.npy')
    flat = pyfits.open('/scratch/dw1519/galex/data/cal/NUV_flat.fits')[0].data
    #factor = 4
    #flat = block_reduce(flat, factor, np.mean)
    frac = cut/full
    frac[np.isnan(frac)] = 0.
    flat = frac*flat
    print flat.shape
    scans = ['0005', '0023', '0050', '0446', '0464', '0482', '0806']
    scans = ['0941', '0995', '2120', '2633','3047']
    #scans = ['0005']
    scans = ['0158', '0689','1247','1508', '1778', '3326', '3587']
    scans = ['3326', '3587']
    scans = ['0212','0302','0500','0743','2291','2453', '2759', '2885', '2930', '3155', '3227', '3497']
    scans = ['0005', '0023', '0050', '0446', '0464', '0482',
             '0806', '0158', '0689','1247','1508', '1778', 
             '3326','0941', '0995', '2120', '2633','3047']
    date = '08-21-2017'#'08-21-2017'
    #model = np.zeros((size,size))
    #data = np.zeros((size,size))
    for scan in scans:
        model = np.zeros((size,size))
        data = np.zeros((size,size))
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
        cut_path = '/scratch/dw1519/galex/co/remove/'

        step=0.01
        asp_solution_list = []
        ix_list = []
        dead_list = []
        cut_list = []
        photon_list = []
        exp_list = []
        for name in name_list:
            photon_list.append(np.load('/scratch/dw1519/galex/data/star_photon/'+name+'_photon.npy'))
            co_data = pyfits.open('{0}/{1}{2}-asprta.fits'.format(asprta, name, asprta_suffix))[1].data
            n = resolution/step
            T = co_data['T']
            ra = co_data['ra']
            dec = co_data['dec']
            roll = co_data['roll']
            t_new = np.arange((T.shape[0]-1)*n)*step+T[0]


            spl_ra = splrep(T, ra)
            spl_dec = splrep(T, dec)
            spl_roll = splrep(T, roll)
            ra_new = splev(t_new, spl_ra)
            dec_new = splev(t_new, spl_dec)
            roll_new = splev(t_new, spl_roll)

            asp_solution_tmp = np.array([t_new, ra_new, dec_new, roll_new]).T

            scst_tmp = pyfits.open('{0}/{1}-scst.fits'.format(scst, name))[1].data
            scst_time = scst_tmp['pktime']
            ix_tmp = np.digitize(asp_solution_tmp[:,0], scst_time)-1
            ix_mask = (ix_tmp>=0) & (ix_tmp<scst_time.shape[0])
            cut = np.load(cut_path+name+'.npy')
            ix_cut = np.digitize(asp_solution_tmp[:,0], cut[:,0]/1000.)-1   
            cut_mask = (ix_cut>=0) & (ix_cut<cut.shape[0])
            id_mask = ix_mask&cut_mask
            ix_tmp = ix_tmp[id_mask]
            ix_cut = ix_cut[id_mask]
            asp_solution_tmp = asp_solution_tmp[id_mask]
            hv = scst_tmp['hvnom_nuv']
            mask = hv[ix_tmp]>0
            ix_tmp = ix_tmp[mask]
            ix_cut = ix_cut[mask]
            asp_solution_tmp = asp_solution_tmp[mask]
            asp_solution_list.append(asp_solution_tmp)
    
            #limit = scst_tmp['t_dead_nuv'].shape[0]-1
            #ix_tmp[ix_tmp>limit] = limit
            dead = scst_tmp['t_dead_nuv']
            tec = scst_tmp['NDCTEC']
            fec = scst_tmp['NDCFEC']
            if name == 'AIS_GAL_SCAN_03497_0002':
                hold = 220000
            else:
                hold = 150000
            fec_mask = fec>hold
            width=1000
            while True:
                ratio_mask = (fec>(hold-width)) & (fec<(hold+width))
                if np.sum(ratio_mask)>10:
                    break
                else:
                    width+=1000
            #ratio = np.mean((1.-dead[ratio_mask])/(tec[ratio_mask]/fec[ratio_mask]))
            ratio = np.median(1.-dead[ratio_mask])/np.median(tec[ratio_mask]/fec[ratio_mask])
            print 'ratio:{0}'.format(ratio)
            dead[fec_mask] = 1.-ratio*tec[fec_mask]/fec[fec_mask]
            dead_tmp = dead[ix_tmp]
            dead_list.append(dead_tmp)
            cut_list.append(cut[ix_cut,1])

        asp = np.concatenate(asp_solution_list, axis=0)
        dead = np.concatenate(dead_list, axis=0)
        cut = np.concatenate(cut_list, axis=0)
        photon_list = np.concatenate(photon_list, axis=0)
        #qmask = photon_list[:,-1]>5
        #photon_list = photon_list[qmask]

        star_list = np.load('/scratch/dw1519/galex/data/star_photon/'+name+'_star.npy')
        print star_list.shape
        stars, ct = np.unique(photon_list[:,0], return_counts=True)
        nuv = star_list[stars.astype(int),-1]
        print stars.shape
        star_mask = (ct>10) & (nuv<18) & (nuv>14.5)
        print np.sum(star_mask)
        stars = stars[star_mask]
        for star in stars:
            star_mask = (photon_list[:,0]==star)
            p = photon_list[star_mask]
            ix = np.digitize(p[:,1]/1000., asp[:,0])-1
            xi, eta = gn.gnomfwd_simple(c_radec[star].ra.deg, c_radec[star].dec.deg, 
                                            asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
            pos = np.array([xi,eta]).T
            pos=((pos/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*size

            Hc, xedge, yedge = np.histogram2d(pos[:,0], pos[:,1], 
                                bins=[size,size], range=([ [0,size],[0,size] ]))


            xi, eta = gn.gnomfwd_simple(c_radec[star].ra.deg, c_radec[star].dec.deg, 
                                            asp[:,1], asp[:,2], -asp[:,3],1/36000.,0.)
            pos = np.array([xi,eta]).T
            pos=((pos/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*size

            H, xedge, yedge = np.histogram2d(pos[:,0], pos[:,1],
                                bins=[size,size], range=([ [0,size],[0,size] ]), weights=step*(1-dead)*cut)

            #A = np.sum(Hc)/np.sum(np.multiply(H,flat))
            #model += A*np.multiply(H,flat)
            data += Hc

            mask = H>0
            if np.sum(mask)>1:
                print np.sum(mask)
                x,y = np.where(mask)
                count = Hc[mask]
                exp = H[mask]
                R = flat[x,y]
                exp_list.append(np.column_stack([star*np.ones(x.shape[0]), x, y, exp, count, R]))
        #np.save('/scratch/dw1519/galex/data/star_photon/'+scan+'_model.npy', model)
        exp_list = np.concatenate(exp_list, axis=0)
        np.save('/scratch/dw1519/galex/data/star_photon/new/'+scan+'_data.npy', data)
        np.save('/scratch/dw1519/galex/data/star_photon/new/'+scan+'_exp.npy', exp_list)

        '''
            mask = H>0
            if np.sum(mask)>1:
                print np.sum(mask)
                x,y = np.where(mask)
                count = Hc[mask]
                exp = H[mask]
                R = flat[x,y]

                print(np.min(exp))
                print(np.max(exp))
                print(np.min(count))
                print(np.max(count))
                print(np.sum(count))

                exp_list.append(np.column_stack([star*np.ones(x.shape[0]), x, y, exp, count, R]))

            inm = count/exp
            plt.plot(y, count, '.k')
            plt.show()
            plt.imshow(Hc)
            plt.show()
            plt.imshow(H)
            plt.show()
            plt.imshow(Hc/H)
            plt.show()
            plt.plot(x, inm/(np.sum(inm)/np.sum(R)), '.k')
            plt.plot(x, R, '.r')
            plt.show()

        exp_list = np.concatenate(exp_list, axis=0)
        np.save('/scratch/dw1519/galex/data/star_photon/'+scan+'_exp.npy', exp_list)
        '''
