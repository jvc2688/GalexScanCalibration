import numpy as np
import re
import glob
import gnomonic as gn
import imagetools
import gc
from astropy.io import fits as pyfits
from scipy import interpolate
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import c3
from matplotlib.colors import LogNorm
import sys


def centroid(t,H,wcs):
    print 'threshold:{0}'.format(t) 

    data = H.byteswap(True).newbyteorder()
    data = data.copy(order='C')
    data = data.byteswap(True).newbyteorder()
    c_data = c3.find_centroid(data, t)
    if c_data is not None:
        cy, cx, max_value, flux = c_data
        cy+=0.5
        cx+=0.5
        centroid = wcs.wcs_pix2world(wcs.sip_foc2pix([[cx, cy]],1),1)[0]
        print(cx,cy)
        print(centroid)
    else:
        centroid = [0.,0.]
        max_value = 0
        flux = 500
    if centroid[0]>1:
        centroid[0] = centroid[0]-360.
    return centroid

def plot_co(co, pixsz, skypos, skyrange, weights=None):
    bound = skyrange[0]/2.
    imsz = imagetools.deg2pix(skypos, skyrange, pixsz).astype(int)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
    foc = wcs.sip_pix2foc(wcs.wcs_world2pix(co[:,0:2],1),1)
    if weights is not None:
        H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]), weights=weights)
    else:
        H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
    h = H.copy()
    #np.flip(h,axis=0)
    '''
    hist_size = [100,100]
    H,xedges,yedges=np.histogram2d(co[:,1], co[:,0],\
                             bins=hist_size, range=([ [-bound,bound],[-bound,bound] ]))
    h = H.copy()
    '''
    plt.imshow(h+10**-10,interpolation='None', cmap=plt.get_cmap('Greys'),\
                 aspect='equal', extent=[bound*3600, -bound*3600, -bound*3600, bound*3600], origin='lower',vmin=0)
                 #norm=LogNorm(10**0,np.max(h)))

    plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
    plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
    plt.xlabel(r'$\Delta gl$')
    plt.ylabel(r'$\Delta gb$')
    plt.xlim(-bound*3600,bound*3600)
    plt.ylim(-bound*3600,bound*3600)

    '''
    data = H.byteswap(True).newbyteorder()
    data = data.copy(order='C')
    data = data.byteswap(True).newbyteorder()
    t=50
    c_data = c3.find_centroid(data, t)
    if c_data is not None:
        cy, cx, max_value, flux = c_data
        c = np.array([-bound+cx*bound*2/hist_size[0], -bound+cy*bound*2/hist_size[0]])
    print c*3600
    '''
    c = centroid(5,H,wcs)
    print(c*3600)
    plt.plot(c[0]*3600,c[1]*3600,'+r', markersize=8, mew=1)

    #plt.show()

def co_rel(data, catalog):
    mask = data[:,-1]>0
    photon = data[mask]
    star = np.unique(photon[:,-1]).astype(int)
    co = np.array([[],[],[],[],[],[]]).T
    print(star.shape)
    for idx in star:
        print(idx)
        star_mask = photon[:,-1]==idx
        print(np.sum(star_mask))
        photon_tmp = photon[star_mask]
        photon_tmp[:,0:2] -= catalog[idx]
        co = np.concatenate([co, photon_tmp])
    return co


def detector2gon(pos_list, size):
  pos_array = np.array(pos_list, dtype='float64')
  pos = (pos_array/50.-0.5)*36000*800*0.001666
  return pos

def get_key(path):
	return int(re.split('_|/|\.', path)[-2])

def dis_correct(data, dis_map, asp, slope, inter):
    ya_mask = data[:,2]>=2
    q_mask = np.ones(data.shape[0], dtype=bool)
    #q_mask = data[:,3]>5
    data = data[ya_mask&q_mask]
    print data.shape
    ix = np.digitize(data[:,0]/1000., asp[:,0])-1

    low_mask = data[:,1]<16
    high_mask = data[:,1]>=16

    coo = data[low_mask, 4:6]
    data[low_mask, 4] -= dis_map['xi_l'](coo)
    data[low_mask, 5] -= dis_map['eta_l'](coo)

    coo = data[high_mask, 4:6]
    data[high_mask, 4] -= dis_map['xi_h'](coo)
    data[high_mask, 5] -= dis_map['eta_h'](coo)

    if slope is not None:
        data[:,5] -= (slope*data[:,2]+inter)*36000*800*0.001666/2400.

    ra, dec = gn.gnomrev_simple(data[:,4], data[:,5], asp[ix,1], asp[ix,2], -asp[ix,3],1/36000.,0.)
    return np.concatenate([np.array([ra,dec]).T, data[:,0:6], asp[ix,1:4]], axis=1)

def load_dis(dis_l, dis_h):
    detector_size = [50,50]

    #dis_map_l = np.load('../data/distortion/194-437-xa_low-centroids.npy')
    #dis_map_h = np.load('../data/distortion/194-437-xa_high-centroids.npy')
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

    return dis_map

def load_ya(ya):
    ya_info = np.load(ya)
    slope = ya_info[0]
    inter = ya_info[1]
    return ya_info

def load_data(name, dis_l, dis_h, ya, c_file, asprta, resolution):
    #name = 'AIS_GAL_SCAN_00023_0001'
    suffix = '-cal-sec'
    #dis_l = '/scratch/dw1519/galex/data/distortion/5-185-xa_low-centroids.npy'
    #dis_h = '/scratch/dw1519/galex/data/distortion/5-185-xa_high-centroids.npy'
    #ya = '/scratch/dw1519/galex/co/co23_1-10/test/ya_fit.npy'

    hdulist = pyfits.open(asprta)
    co_data = hdulist[1].data
    n=resolution/0.005 #100
    T = co_data['T']
    ra = co_data['ra']
    dec = co_data['dec']
    roll = co_data['roll']
    t_new = np.arange((T.shape[0]-1)*n)*0.005+T[0]
    ra_new = np.interp(t_new, T, ra)
    dec_new = np.interp(t_new, T, dec)
    roll_new = np.interp(t_new, T, roll)
    asp = np.array([t_new, ra_new, dec_new, roll_new]).T
    #np.save('../data/photon_list/%s_asp_fix.npy'%name, asp_cal)
    hdulist.close()
    
    #asp = np.load('../data/photon_list/{0}{1}-correct_asp_new.npy'.format(name, suffix))
    scst_file = pyfits.open('../AIS_GAL_SCAN/scst/{0}-scst.fits'.format(name))
    scst_data = scst_file[1].data
    scst_time = scst_data['pktime'].copy()
    hv = scst_data['hvnom_nuv'].copy()
    scst_file.close()

    dis_map = load_dis(dis_l, dis_h)
    if ya is not None:
        ya_info = load_ya(ya)
    else:
        ya_info = [None, None]


    tycho2 = pyfits.open('../data/tycho2.fits')[1].data
    df = load_catalog(c_file)
    c = SkyCoord(df['gl']*u.degree, df['gb']*u.degree, frame='galactic')
    catalog = SkyCoord(tycho2['Glon']*u.degree, tycho2['Glat']*u.degree, frame='galactic')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
    mask=d2d<0.001*u.degree
    photon_list =  sorted(glob.glob("../data/"+name+suffix+"/split/*.csv"), key=get_key)
    data_list = []
    for f in photon_list:
        print(f)
        data = np.loadtxt(f,delimiter=',', usecols=[0,3,4,5,6,7])
        if len(data.shape)<2:
            continue
        scst_ix = np.digitize(data[:,0]/1000., scst_time)
        ix_mask = scst_ix<scst_time.shape[0]
        scst_ix = scst_ix[ix_mask]
        data = data[ix_mask]
        print(data[-1,0]/1000.)
        print(scst_time[scst_ix[-1]])
        print(hv[scst_ix[-1]])
        data = data[hv[scst_ix]>0]
        if len(data)==0:
            print('skip')
            continue
        data = dis_correct(data, dis_map, asp, ya_info[0], ya_info[1])
        sky_data = SkyCoord(data[:,0:2], unit='deg', frame='fk5', equinox='J2000.0')
        gal = sky_data.transform_to('galactic')
        idxc, idxcatalog, d2d_p, d3d_p = gal.search_around_sky(catalog[idx[mask]], 0.0045*u.deg)
        catalog_id = np.zeros((gal.shape[0], 3))-1
        #catalog_id[idxcatalog,0] = 
        if len(idxc)>0:
            xi, eta = gn.gnomfwd_simple(tycho2['RAJ2000'][idx[mask][idxc]], tycho2['DEJ2000'][idx[mask][idxc]],
                     data[idxcatalog,-3], data[idxcatalog,-2], -data[idxcatalog,-1], 1/36000., 0)
            #print xi.shape
            #print eta.shape
            #print data[idxcatalog,6].shape
            #print (data[idxcatalog,7]-eta).shape        
            catalog_id[idxcatalog] = np.column_stack([idx[mask][idxc], (data[idxcatalog,6]-xi),\
                                    (data[idxcatalog,7]-eta)])
        else:
            print 'no star'
        #catalog_id = np.array([catalog_id]).T
        data = np.concatenate([np.array([gal.l.deg]).T, np.array([gal.b.deg]).T, data[:,2:8],catalog_id], axis=1)
        data_list.append(data)
    gc.collect()

    return np.concatenate(data_list, axis=0)

def load_data_one(i):
    name = 'AIS_GAL_SCAN_00023_0001'
    suffix = '-cal-sec'
    dis_l = '/scratch/dw1519/galex/data/distortion/5-185-xa_low-centroids.npy'
    dis_h = '/scratch/dw1519/galex/data/distortion/5-185-xa_high-centroids.npy'
    ya = '/scratch/dw1519/galex/co/co23_1-10/test/ya_fit.npy'
    
    asp = np.load('../data/photon_list/{0}{1}-correct_asp_new.npy'.format(name, suffix))
    scst_file = pyfits.open('../AIS_GAL_SCAN/scst/{0}-scst.fits'.format(name))
    scst_data = scst_file[1].data
    scst_time = scst_data['pktime'].copy()
    hv = scst_data['hvnom_nuv'].copy()
    scst_file.close()

    dis_map = load_dis(dis_l, dis_h)
    ya_info = load_ya(ya)
    photon_list =  sorted(glob.glob("../data/"+name+suffix+"/split/*.csv"), key=get_key)
    data_list = []
    f = photon_list[i]
    data = np.loadtxt(f,delimiter=',', usecols=[0,3,4,5,6,7])
    scst_ix = np.digitize(data[:,0]/1000., scst_time)
    ix_mask = scst_ix<scst_time.shape[0]
    scst_ix = scst_ix[ix_mask]
    data = data[ix_mask]
    print(data[-1,0]/1000.)
    print(scst_time[scst_ix[-1]])
    print(hv[scst_ix[-1]])
    data = data[hv[scst_ix]>0]
    if len(data)==0:
        print('skip')
    data = dis_correct(data, dis_map, asp, ya_info[0], ya_info[1])

    return data

def load_catalog(c_file):
    column = ["NUMBER","X_IMAGE","Y_IMAGE","ALPHA_J2000","DELTA_J2000",\
            "FLUX_AUTO","FLUXERR_AUTO","FLUX_APER","A_IMAGE","B_IMAGE",
            "THETA_IMAGE","FWHM_IMAGE","x_new","y_new","nuv","gl","gb"]
    df = pd.read_csv(c_file,skiprows=[0,1,2,3],sep='\s+', names=column)
    #df = pd.read_csv('../data/extraction/Starcat-07-22-2017/starcat_0023mapweight_fec_fwhm.txt',\
    # skiprows=[0,1,2,3],sep='\s+', names=column)
    return df


if __name__ == '__main__':
    if False:
        name = sys.argv[1]
        dis_l = sys.argv[2]
        dis_h = sys.argv[3]
        ya = sys.argv[4]
        c_file = sys.argv[5]
        asprta = sys.argv[6]
        resolution = float(sys.argv[7])

        data = load_data(name, dis_l, dis_h, ya, c_file, asprta, resolution)
        np.save('../data/'+name+'-cal-sec/photon_match_d_072217.npy', data)

        '''
        tycho2 = pyfits.open('../data/tycho2.fits')[1].data
        catalog = np.array([tycho2['Glon'], tycho2['Glat']]).T
        co = co_rel(data, catalog)
        np.save('../data/'+name+'-cal-sec/co_match_072217.npy', co)
        '''

    if True:
        name = sys.argv[1]
        dis_l = sys.argv[2]
        dis_h = sys.argv[3]
        ya = None#sys.argv[4]
        c_file = sys.argv[5]
        asprta = sys.argv[6]
        resolution = float(sys.argv[7])

        data = load_data(name, dis_l, dis_h, ya, c_file, asprta, resolution)
        np.save('../data/'+name+'-cal-sec/photon_match_d_081017.npy', data)

        '''
        tycho2 = pyfits.open('../data/tycho2.fits')[1].data
        catalog = np.array([tycho2['Glon'], tycho2['Glat']]).T
        co = co_rel(data, catalog)
        np.save('../data/'+name+'-cal-sec/co_match_072217.npy', co)
        '''