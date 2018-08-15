import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits
import imagetools
import gnomonic as gn
import pos_range
import scipy
import math
import threading
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from multiprocessing import Process
from multiprocessing import Manager
import sys
import os
import gc
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
import argparse
from astropy import units as u



def exp_dead_new(file_num, name_file, imsz, wcs, hrflat, foc_list, asp_solution, dead, cut, step, return_dict):
    print imsz
    count = np.zeros(imsz)

    x_lim = imsz[0]
    y_lim = imsz[1]

    length = hrflat.shape[0]
    half_len = length/2.
    print half_len
    l = imsz[0]/10
    start = foc_list[0,1]-half_len
    print foc_list.shape
    print start.shape
    ox = np.repeat(np.arange(l)+start,imsz[1])
    oy = np.tile(np.arange(imsz[1]),l)
    omask = ox<imsz[0]
    ox = ox[omask]
    oy = oy[omask]
    gl,gb = wcs.all_pix2world(oy,ox,0)
    c = SkyCoord(gl*u.degree, gb*u.degree, frame='galactic')
    rd = c.transform_to(FK5)
    for i in range(asp_solution.shape[0]):
        foc = foc_list[i,:]#wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
        if (foc[1]+half_len)>=(start+l):
            print 'update'
            start = foc[1]-half_len
            ox = np.repeat(np.arange(l)+start,imsz[1])
            oy = np.tile(np.arange(imsz[1]),l)
            omask = ox<imsz[0]
            ox = ox[omask]
            oy = oy[omask]
            gl,gb = wcs.all_pix2world(oy,ox,0)
            c = SkyCoord(gl*u.degree, gb*u.degree, frame='galactic')
            rd = c.transform_to(FK5)
        fmask = (ox>=(foc[1]-length/2)) & (ox<(foc[1]+length/2)) & (oy>=(foc[0]-length/2)) & (oy<(foc[0]+length/2))
        x = ox[fmask]
        y = oy[fmask]
        xi, eta = gn.gnomfwd_simple(rd.ra.deg[fmask], rd.dec.deg[fmask], 
                                        asp_solution[i,1], asp_solution[i,2], -asp_solution[i,3],1/36000.,0.)
        px = ((xi/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*length
        py = ((eta/36000.)/(1.25/2.)*(1.25/(800* 0.001666))+1.)/2.*length
        pmask = (px>=0) & (px<length) & (py>=0) & (py<length)

        count[x[pmask].astype(int),y[pmask].astype(int)] += \
            hrflat[px[pmask].astype(int),py[pmask].astype(int)]*step*(1-dead[i])*cut[i]
        if i%100==0:
            with open('/scratch/dw1519/galex/fits/scan_map/%s_gal_sec_exp_tmp%d.dat'%(name_file, file_num),'w') as f:
                f.write('%d'%i)
            print i
    print '%d done'%file_num
    return_dict[file_num] = count
    #np.save('../fits/scan_map/%s_gal_sec_exp_tmp%d.npy'%(name_file, file_num), count)



if __name__ == '__main__':

#process
    if True:
        parser = argparse.ArgumentParser(description='galex scan count map')
        parser.add_argument('scan_name', nargs=1, help="scan name")
        parser.add_argument('np', nargs=1, type=int, help="number of process")
        parser.add_argument('guide_path', nargs=1, help="path to the guide asp file")
        parser.add_argument('guide_suffix', nargs=1, help="suffix of the guide asp file")
        parser.add_argument('asprta', nargs=1, help="path to the asprta file")
        parser.add_argument('asprta_suffix', nargs=1, help="suffix the asprta file")
        parser.add_argument('asprta_resolution', nargs=1, type=float, help="resolution of asprta")
        parser.add_argument('scst', nargs=1, help="path to the scst file")
        parser.add_argument('out_path', nargs=1, help="path to the output directory")
        parser.add_argument('-c', '--cut', nargs=1, help="path to cut correction file")

        args = parser.parse_args()

        name_file = args.scan_name[0]

        guide_path = args.guide_path[0]
        guide_suffix = args.guide_suffix[0].lstrip()

        asprta_path = args.asprta[0]
        asprta_suffix = args.asprta_suffix[0].lstrip()
        resolution = args.asprta_resolution[0]
        scst = args.scst[0]

        num_t = args.np[0]

        out_path = args.out_path[0]

        if args.cut is not None:
            cut_path = args.cut[0]
        else:
            cut_path = None

        print("scan name: {0}".format(name_file))
        print("number of process: {0}".format(num_t))
        print("guide asp path: {0}".format(guide_path))
        print("guide asp suffix: {0}".format(guide_suffix))
        print("asprta: {0}".format(asprta_path))
        print("asprta suffix: {0}".format(asprta_suffix))
        print("resolution: {0}".format(resolution))
        print("scst: {0}".format(scst))
        print("output: {0}".format(out_path))
        print("cut: {0}".format(cut_path))
        
        with open('../name_scan/%s'%name_file) as f:
            name_list = f.read().splitlines()
        print name_list

        pixsz = 0.0005555555556 #0.000416666666666667 #0.00166667
        step=0.05
        asp_solution_list = []
        ix_list = []
        dead_list = []
        cut_list = []
        corr_list = []
        gal_l = []
        for name in name_list:
            #try:
            #    asp_solution_tmp = np.load('../data/photon_list/{0}{1}-correct_asp_new_bs.npy'.format(name, suffix))
            #except IOError:
            co_data = pyfits.open('{0}/{1}{2}-asprta.fits'.format(asprta_path, name, asprta_suffix))[1].data
            T = co_data['T']
            ra = co_data['ra']
            dec = co_data['dec']
            roll = co_data['roll']
            n = resolution/step
            t_new = np.arange((T.shape[0]-1)*n)*step+T[0]

            '''
            f_ra = interp1d(T, ra, kind='cubic')
            f_dec = interp1d(T, dec, kind='cubic')
            f_roll = interp1d(T, roll, kind='cubic')

            ra_new = f_ra(t_new)
            dec_new = f_dec(t_new)
            roll_new = f_roll(t_new)
            '''
            
            spl_ra = splrep(T, ra)
            spl_dec = splrep(T, dec)
            spl_roll = splrep(T, roll)
            ra_new = splev(t_new, spl_ra)
            dec_new = splev(t_new, spl_dec)
            roll_new = splev(t_new, spl_roll)
        
            '''
            ra_new = np.interp(t_new, T, ra)
            dec_new = np.interp(t_new, T, dec)
            roll_new = np.interp(t_new, T, roll)
            '''
            asp_solution_tmp = np.array([t_new, ra_new, dec_new, roll_new]).T
            #    np.save('../data/photon_list/{0}{1}-correct_asp_new_bs.npy'.format(name, suffix), asp_solution_tmp)
            #ix_tmp = (np.round(asp_solution_tmp[:,0]-asp_solution_tmp[0,0])+1).astype(int)
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

            #asp_old = np.load('../data/photon_list/%s_asp.npy'%name)
            #sky_data = SkyCoord(asp_old[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
            #calculate the image size and location
            #asp_map = np.load('{0}/{1}{2}_asp.npy'.format(guide_path, name, guide_suffix))
            asp_map = np.load('{0}/{1}_asp.npy'.format(guide_path, name))
            sky_data = SkyCoord(asp_map[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
            gal = sky_data.transform_to(Galactic)
            asprta = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
            gal_l.append(np.mean(asprta[:,0]))

        asp_solution = np.concatenate(asp_solution_list, axis=0)
        sky_data = SkyCoord(asp_solution[:,1:3], unit='deg', frame=FK5, equinox='J2000.0')
        gal = sky_data.transform_to(Galactic)
        coo = np.concatenate((np.array([gal.l.deg]).T, np.array([gal.b.deg]).T), axis=1)
        
        dead = np.concatenate(dead_list, axis=0)
        cut = np.concatenate(cut_list, axis=0)
        
        flat_hdu = pyfits.open('../data/cal/NUV_flat.fits')
        flat =flat_hdu[0].data
        size=flat_hdu[0].header['CDELT2']
        flat = np.load('/scratch/dw1519/galex/data/star_photon/super_iterate3/flat.npy')
        resample=3 #4
        hrflat = scipy.ndimage.interpolation.zoom(flat, resample, order=0,prefilter=False)#scipy.ndimage.interpolation.zoom(flat,size/pixsz,order=1,prefilter=False)
        print hrflat.shape
        #hrflat = np.swapaxes(hrflat,0,1)
        #print hrflat.shape
        #hrflat = scipy.ndimage.interpolation.rotate(hrflat,-23.4,reshape=False,order=1,prefilter=False)
        #print hrflat.shape
        skypos = [np.mean(gal_l), 0.]
        #skyrange = [np.max(gal_l)-np.min(gal_l)+1.5, 21.5]
        skyrange = [1.55, 21.5]
        tranges = []
        trange = [0, 0]
        imsz = imagetools.deg2pix_g(skypos, skyrange, pixsz).astype(int)
        wcs = imagetools.define_wcs_g(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)

        focs = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,0),0)
        idx = np.argsort(focs[:,1])

        focs = focs[idx]
        asp_solution = asp_solution[idx]
        dead = dead[idx]
        cut = cut[idx]

        foc_list = np.array_split(focs, num_t)
        asp_solutions = np.array_split(asp_solution, num_t)
        deads = np.array_split(dead, num_t)
        cuts = np.array_split(cut, num_t)


        np.save('/scratch/dw1519/galex/fits/scan_map/dead/dead_{1}.npy'.format(out_path, name_file), dead)
        plt.plot(dead,'.')
        plt.savefig('/scratch/dw1519/galex/fits/scan_map/dead/dead_{1}.pdf'.format(out_path, name_file))
        plt.clf()


        manager = Manager()
        return_dict = manager.dict()
        p_list = []

        for file_num in range(0, num_t):
            p = Process(target=exp_dead_new, args=(file_num, name_file, imsz, wcs, hrflat, foc_list[file_num], asp_solutions[file_num], deads[file_num], cuts[file_num], step, return_dict))
            p.start()
            p_list.append(p)

        for p in p_list:
            p.join()
        print 'all done'
        
        count = np.zeros(imsz)

        for i in range(0, num_t):
            #count += np.load('{0}/{1}_gal_sec_exp_tmp{2}.npy'.format(out_path, name_file, i)) #return_dict[i]
            count += return_dict[i]

        tranges.append(trange)
        hdu = pyfits.PrimaryHDU(count)
        hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('{0}/count_map_{1}_exp.fits'.format(out_path, name_file), clobber=False)


        for i in range(num_t):
            #os.remove('{0}/{1}_gal_sec_exp_tmp{2}.npy'.format(out_path, name_file, i))
            os.remove('/scratch/dw1519/galex/fits/scan_map/{1}_gal_sec_exp_tmp{2}.dat'.format(out_path, name_file, i))



