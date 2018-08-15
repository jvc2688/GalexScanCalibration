import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import c3
from astropy.io import fits as pyfits
import numpy as np
from astropy import wcs as pywcs
import sys
import aplpy
import os
from sklearn.neighbors import KernelDensity
import csv
import math
import asp_cal
import glob
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.io.fits import update
import re
from scipy.interpolate import splev, splrep

def _find_centroid(filename):
    try:
        hdulist = pyfits.open(filename)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
        hdulist = None
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise
    if hdulist is not None:
        w = pywcs.WCS(hdulist[0].header, hdulist)
        data = hdulist[0].data
        data = data.byteswap(True).newbyteorder()
        cy, cx = c3.find_centroid(data)
    else:
        return None
    return cx,cy

def find_centroid(filename):
    try:
        hdulist = pyfits.open(filename)
    except IOError as e:
        print "I/O error({0}): {1}: {2}".format(e.errno, e.strerror, filename)
        hdulist = None
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise
    if hdulist is not None:
        hdulist = pyfits.open(filename)
        w = pywcs.WCS(hdulist[0].header, hdulist)
        data = hdulist[0].data
        data = data.byteswap(True).newbyteorder()
        cy, cx = c3.find_centroid(data)
        centroid = w.wcs_pix2world(w.sip_foc2pix([[cx, cy]],1),1)[0]
        if centroid[0]>1:
            centroid[0] = centroid[0]-360.
    else:
        centroid = [0,0]
    return centroid

def get_centers(initial, final):
    centers = []
    for i in range(initial, final+1):
        filename = '../fits/co/right/co_map%d_%d_zoom_nocal.fits'%(i,i+1)
        centroid = find_centroid(filename)
        centers.append(centroid)
    centers = np.array(centers)
    return centers

def get_centers_half(initial, final):
    centers = []
    for i in range(initial, final+1):
        for j in range(10):
            filename = '../fits/test/co_map%d_%d_10.fits'%(i,j)
            centroid = find_centroid(filename)
            centers.append(centroid)
    centers = np.array(centers)
    return centers

def corr_plot(centroid, filename, title):

    print filename
    #centroid = find_centroid(filename)
    #centroid = np.load('../data/offsets300-899_half_r.npy')[7]
    print centroid
    fig = aplpy.FITSFigure(filename)
    fig.add_label(centroid[0], centroid[1], 'X', color='red')
    fig.show_grayscale(invert=True)
    fig.tick_labels.set_xformat('d.ddd')
    fig.tick_labels.set_yformat('d.ddd')
    fig.recenter(0., 0., radius=0.01)
    fig.add_grid()
    fig.set_title(title)
    basename = os.path.basename(filename)
    preffix, ext = os.path.splitext(basename)
    fig.save('../plots/corr_10/%s.png'%preffix)

def load_data(filename):
    data = []
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if math.fabs(float(row[0]))>0.008 or math.fabs(float(row[1]))>0.008:
                pass
            else:
                data.append(row)
        csvfile.close()
    return data

def moving_stat(data, out_mask, half_win=10):
    moving_mask = np.zeros(data.shape)
    moving_mask[0:2*half_win+1] = 1
    mean = np.zeros(data.shape)
    median = np.zeros(data.shape)
    abs_dev = np.zeros(data.shape)
    std = np.zeros(data.shape)
    z = np.zeros(data.shape)
    mz = np.zeros(data.shape)

    for i in range(half_win):
        if out_mask[i] == 1:
            std[i] = 1
            abs_dev[i] = 1
        else:
            tmp_out_mask = -out_mask[:half_win+i+1]+1
            #print i, data[:half_win+i+1][tmp_out_mask>0].shape
            mean[i] = np.mean(data[:half_win+i+1][tmp_out_mask>0], axis=0)
            std[i] = np.std(data[:half_win+i+1][tmp_out_mask>0], axis=0)

            median[i] = np.median(data[:half_win+i+1][tmp_out_mask>0], axis=0)
            abs_dev[i] = np.median(np.absolute(data[:half_win+i+1][tmp_out_mask>0]-median[i]), axis=0)

        if out_mask[-i-1] == 1:
            std[-i-1] = 1
            abs_dev[-i-1] =1
        else:
            tmp_out_mask = -out_mask[-half_win-i-1:]+1
            median[-i-1] = np.median(data[-half_win-i-1:][tmp_out_mask>0], axis=0)
            abs_dev[-i-1] = np.median(np.absolute(data[-half_win-i-1:][tmp_out_mask>0]-median[-i-1]), axis=0)

            mean[-i-1] = np.mean(data[-half_win-i-1:][tmp_out_mask>0], axis=0)
            std[-i-1] = np.std(data[-half_win-i-1:][tmp_out_mask>0], axis=0)
            #print -i-1, data[-half_win-i-1:][tmp_out_mask>0].shape

    for i in range(data.shape[0]-2*half_win):
        if out_mask[half_win+i] == 1:
            std[half_win+i] = 1
            abs_dev[half_win+i] =1
        moving_mask_tmp = np.roll(moving_mask, i)
        tmp_out_mask = -out_mask[moving_mask_tmp>0]+1
        mean[half_win+i] = np.mean(data[moving_mask_tmp>0][tmp_out_mask>0], axis=0)
        std[half_win+i] = np.std(data[moving_mask_tmp>0][tmp_out_mask>0], axis=0)
        median[half_win+i] = np.median(data[moving_mask_tmp>0][tmp_out_mask>0], axis=0)
        abs_dev[half_win+i] = np.median(np.absolute(data[moving_mask_tmp>0][tmp_out_mask>0]-median[half_win+i]), axis=0)

        #print half_win+i, data[moving_mask_tmp>0][tmp_out_mask>0].shape

    z = np.absolute((data - mean)/std)
    mz = np.absolute(0.6745*(data-median)/abs_dev)

    return z, mz

def generate_zero_offsets(name):
    print name
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    initial = 1
    final = hdulist[1].data['T'].shape[0]-1

    centers = []
    for i in range(initial, final+1):
        centers.append(np.load('../data/%s/cata/centroids_photon%d.npy'%(name, i)))
    centroids = np.concatenate(centers, axis=0)

    print centroids.shape

    np.save('../data/%s/cata/offsets%d_10_new_photon.npy'%(name, initial), centroids)

    asp_cal.interpolate_offsets(name)

    output = "../plots/%s/cata/output.csv"%(name)
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
        os.makedirs(dir)

def generate_first_offsets(name):
    print name
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    initial = 1
    final = hdulist[1].data['T'].shape[0]-1

    centers = []
    center_time = []
    for i in range(initial, final+1):
        c = np.load('../data/%s/cata/centroids_rot%d.npy'%(name, i))
        #if c.shape == (1,3):
        #   c = c[:,:2]
        centers.append(c)
        center_time.append(np.load('../data/%s/cata/time_rot%d.npy'%(name, i)))
        print c.shape
    centroids = np.concatenate(centers, axis=0)
    time =  np.concatenate(center_time, axis=0)

    print centroids.shape

    out_mask = np.zeros(centroids.shape[0])

    z, mz = moving_stat(centroids[:,0], out_mask, half_win=100)
    outliers = np.zeros(centroids.shape[0])
    outliers[mz>3.5] = 1
    outliers[out_mask>0] = 1
    outliers = outliers>0
    index = np.arange(centroids.shape[0])
    centroids[outliers, 0] = np.interp(index[outliers], index[~outliers], centroids[~outliers,0])

    z, mz = moving_stat(centroids[:,1], out_mask, half_win=100)
    outliers = np.zeros(centroids.shape[0])
    outliers[mz>3.5] = 1
    outliers[out_mask>0] = 1
    outliers = outliers>0
    index = np.arange(centroids.shape[0])
    centroids[outliers, 1] = np.interp(index[outliers], index[~outliers], centroids[~outliers,1])


    output = "../plots/%s/cata/output.csv"%(name)
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
        os.makedirs(dir)

    plt.plot(centroids[:,0], '.b')
    plt.savefig('../plots/%s/cata/offsets_10_new_half.pdf'%name, dpi=190)
    plt.clf()

    np.save('../data/%s/cata/time%d_10_new_half.npy'%(name, initial), time)
    np.save('../data/%s/cata/offsets%d_10_new_half.npy'%(name, initial), centroids)

    co_data = hdulist[1].data
    T = co_data['T']
    ra = co_data['ra']
    dec = co_data['dec']
    roll = co_data['roll']
    ra_new = np.interp(time, T, ra) - centroids[:,0]
    dec_new = np.interp(time, T, dec) - centroids[:,1]
    roll_new = np.interp(time, T, roll) - centroids[:,2]

    other = np.zeros((time.shape[0], 8))
    array = np.concatenate([np.array([time, ra_new, dec_new, roll_new]).T, other], axis=1)
    data = np.core.records.fromarrays(array.transpose(), dtype=[('T', float), ('RA', float), ('DEC', float), ('ROLL', float),\
                 ('STATUS_FLAG', int), ('ROLL_RAD', float), ('X', float), ('Y', float), ('Z', float), ('XDOT', float), ('YDOT', float), ('ZDOT', float)])

    new_file = '../AIS_GAL_SCAN/asprta/%s-cal-asprta.fits'%(name)
    os.system('cp ../AIS_GAL_SCAN/asprta/%s-asprta.fits ../AIS_GAL_SCAN/asprta/%s-cal-asprta.fits'%(name, name))
    update(new_file, data, 1)
    hdu = pyfits.open(new_file)
    print hdu[1].data['RA'].shape
    print hdu[1].data['DEC'].shape
    hdu.close()

    #asp_cal.interpolate_offsets(name, 1., centroids)
    tmp_files = glob.glob("../data/%s/cata/centroids_rot*"%name)
    for tmp_file in tmp_files:
        os.remove(tmp_file)

    tmp_files = glob.glob("../data/%s/cata/time_rot*"%name)
    for tmp_file in tmp_files:
        os.remove(tmp_file)

    '''
    print centroids.shape
    os.system('cp ../AIS_GAL_SCAN/asprta/%s-asprta.fits ../AIS_GAL_SCAN/asprta/%s-asprta-cal.fits'%(name, name))
    hdu = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta-cal.fits'%(name), 'update')
    hdu[1].data['RA'][1:] -= centroids[:,0]
    hdu[1].data['DEC'][1:] -= centroids[:,1]
    print hdu[1].data['RA'][1:]
    print hdu[1].data['DEC'][1:]
    hdu.flush()
    hdu.close()
    '''

def generate_sec_offsets(name):
    print name
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    initial = 1
    final = hdulist[1].data['T'].shape[0]-1

    centers = []
    for i in range(initial, final+1):
        centers.append(np.load('../data/%s/cata/centroids_sec%d.npy'%(name, i)))
    centroids = np.concatenate(centers, axis=0)

    print centroids.shape
    plt.plot(centroids[:,0], '.b')
    plt.savefig('../plots/%s/cata/offsets_10_new_sec.pdf'%name, dpi=190)
    plt.clf()

    np.save('../data/%s/cata/offsets%d_10_new_sec.npy'%(name, initial), centroids)

    asp_cal.secondary_cal(name)

    tmp_files = glob.glob("../data/%s/cata/centroids_sec*"%name)
    print tmp_files
    for tmp_file in tmp_files:
        os.remove(tmp_file)

def generate_sec_offsets_new(name):
    print name
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    initial = 0
    final = hdulist[1].data['T'].shape[0]-1

    centers = []
    center_time = []
    for i in range(initial, final+1):
        #center = np.load('../data/%s/cata/centroids_rot%d.npy'%(name, i))
        #if center.shape!=(2,3):
        #   print i, center.shape
        centers.append(np.load('../data/%s/cata/centroids_rot%d.npy'%(name, i)))
        center_time.append(np.load('../data/%s/cata/time_rot%d.npy'%(name, i)))
    centroids = np.concatenate(centers, axis=0)
    time =  np.concatenate(center_time, axis=0)

    print centroids.shape
    print time.shape
    output = '../plots/%s/cata/offsets_10_new_sec.pdf'%name
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
      os.makedirs(dir)
    plt.plot(centroids[:,0], '.b')
    plt.savefig('../plots/%s/cata/offsets_10_new_sec.pdf'%name, dpi=190)
    plt.clf()

    np.save('../data/%s/cata/offsets%d_10_new_sec.npy'%(name, initial), centroids)
    np.save('../data/%s/cata/time%d_10_new_sec.npy'%(name, initial), time)
    
    co_data = hdulist[1].data
    T = co_data['T']
    ra = co_data['ra']
    dec = co_data['dec']
    roll = co_data['roll']
    ra_new = np.interp(time, T, ra) - centroids[:,0]
    dec_new = np.interp(time, T, dec) - centroids[:,1]
    roll_new = np.interp(time, T, roll) - centroids[:,2]

    other = np.zeros((time.shape[0], 8))
    array = np.concatenate([np.array([time, ra_new, dec_new, roll_new]).T, other], axis=1)
    data = np.core.records.fromarrays(array.transpose(), dtype=[('T', float), ('RA', float), ('DEC', float), ('ROLL', float),\
                 ('STATUS_FLAG', int), ('ROLL_RAD', float), ('X', float), ('Y', float), ('Z', float), ('XDOT', float), ('YDOT', float), ('ZDOT', float)])

    new_file = '../AIS_GAL_SCAN/asprta/%s-sec-asprta.fits'%(name)
    os.system('cp ../AIS_GAL_SCAN/asprta/%s-asprta.fits ../AIS_GAL_SCAN/asprta/%s-sec-asprta.fits'%(name, name))
    update(new_file, data, 1)
    hdu = pyfits.open(new_file)
    print hdu[1].data['RA'].shape
    print hdu[1].data['DEC'].shape
    hdu.close()

    tmp_files = glob.glob("../data/%s/cata/centroids_rot*"%name)
    for tmp_file in tmp_files:
        os.remove(tmp_file)

    tmp_files = glob.glob("../data/%s/cata/time_rot*"%name)
    for tmp_file in tmp_files:
        os.remove(tmp_file)

def generate_new_offsets_new(name, asprta, suffix, tmp_dir, num_p):
    print name
    #centers = []
    #center_time = []
    '''
    hdulist = pyfits.open('../AIS_GAL_SCAN/asprta/%s-asprta.fits'%name)
    initial = 0
    final = hdulist[1].data['T'].shape[0]-1
    hdulist.close()


    for i in range(initial, final+1):
        #center = np.load('../data/%s/cata/centroids_rot%d.npy'%(name, i))
        #if center.shape!=(2,3):
        #   print i, center.shape
        centers.append(np.load(tmp_dir+'/centroids_rot%d.npy'%(i)))
        center_time.append(np.load(tmp_dir+'/time_rot%d.npy'%(i)))
    '''
    '''
    for i in range(num_p):
        centers.append(np.load(tmp_dir+'/centroids_tmp%d.npy'%(i)))
        center_time.append(np.load(tmp_dir+'/time_tmp%d.npy'%(i)))
    centroids = np.concatenate(centers, axis=0)
    time =  np.concatenate(center_time, axis=0)
    '''
    try:
        centroids = np.load(tmp_dir+'/offsets_%s.npy'%(suffix))
        time = np.load(tmp_dir+'/time_%s.npy'%(suffix))
    except IOError:
        try:
            centroids = np.load(tmp_dir+'/offsets0_%s.npy'%(suffix))
            time = np.load(tmp_dir+'/time0_%s.npy'%(suffix))
        except IOError:
            print 'no file'
            return 0

    ct1 = np.trim_zeros(centroids[:,0],'f')
    lz = centroids[:,0].shape[0]-ct1.shape[0]
    tz = np.trim_zeros(centroids[:,0],'b').shape[0]
    data_mask = np.zeros(centroids.shape[0])
    data_mask[lz:tz] = 1

    lim=10
    out_mask = np.zeros(tz-lz)
    z, mz = moving_stat(centroids[data_mask>0,0],out_mask,half_win=100)
    outlier = mz>lim

    z, mz = moving_stat(centroids[data_mask>0,1],out_mask,half_win=100)
    outlier = outlier | (mz>lim)

    c_new = np.zeros(centroids.shape)
    spl_ra = splrep(time[data_mask>0][~outlier], centroids[data_mask>0,0][~outlier])
    spl_dec = splrep(time[data_mask>0][~outlier], centroids[data_mask>0,1][~outlier])
    c_new[data_mask>0,0] = splev(time[data_mask>0], spl_ra)
    c_new[data_mask>0,1] = splev(time[data_mask>0], spl_dec)

    print centroids.shape
    print time.shape
    output = '../plots/0/%s-%s.pdf'%(name, suffix)
    dir = os.path.dirname(output)
    if not os.path.exists(dir):
      os.makedirs(dir)
    plt.plot(centroids[:,0], '.k')
    plt.plot(c_new[:,0], '.r')
    plt.savefig('../plots/0/%s-%s_ra.pdf'%(name, suffix), dpi=190)
    plt.clf()

    plt.plot(centroids[:,1], '.k')
    plt.plot(c_new[:,1], '.r')
    plt.savefig('../plots/0/%s-%s_dec.pdf'%(name, suffix), dpi=190)
    plt.clf()

    #np.save(tmp_dir+'/offsets_%s.npy'%(suffix), centroids)
    #np.save(tmp_dir+'/time_%s.npy'%(suffix), time)
    
    hdulist = pyfits.open(asprta)
    co_data = hdulist[1].data
    T = co_data['T']
    ra = co_data['ra']
    dec = co_data['dec']
    roll = co_data['roll']
    ra_new = np.interp(time, T, ra) - c_new[:,0]
    dec_new = np.interp(time, T, dec) - c_new[:,1]
    roll_new = np.interp(time, T, roll) - c_new[:,2]

    other = np.zeros((time.shape[0], 8))
    array = np.concatenate([np.array([time, ra_new, dec_new, roll_new]).T, other], axis=1)
    data = np.core.records.fromarrays(array.transpose(), dtype=[('T', float), ('RA', float), ('DEC', float), ('ROLL', float),\
                 ('STATUS_FLAG', int), ('ROLL_RAD', float), ('X', float), ('Y', float), ('Z', float), ('XDOT', float), ('YDOT', float), ('ZDOT', float)])

    new_file = re.split('-asprta.fits', asprta)[0]+'-'+suffix+'-asprta.fits'
    os.system('cp {0} {1}'.format(asprta, new_file))
    update(new_file, data, 1)
    hdu = pyfits.open(new_file)
    print hdu[1].data['RA'].shape
    print hdu[1].data['DEC'].shape
    hdu.close()



    
if __name__ == '__main__':
    if True:
        name = sys.argv[1]
        suffix = sys.argv[2]
        generate_new_offsets_new(name, suffix)
