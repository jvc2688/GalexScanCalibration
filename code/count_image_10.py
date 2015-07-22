import csv
import imagetools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc
from astropy.io import fits as pyfits
import pos_range
import math
import spilt_csv as spi
import astropy.coordinates as pycoo

def cal_photon(offsets, initial, time, row, step):
  if time<initial:
    return np.array(row, dtype='float64')
  else:
    index = int(math.floor((time-initial)*step/1000.-1))
    row = np.array(row, dtype='float64')
    row[-3:-1] += offsets[index]
    return row 

def cal_photon_r(offsets, initial, time, row, step):
  off_len = offsets.shape[0]
  if time<initial:
    return np.array(row, dtype='float64')
  else:
    index = int(math.floor((time-initial)*step/1000.))
    row = np.array(row, dtype='float64')
    if index >= off_len:
      return row
    row[-3:-1] -= offsets[index]
    return row 

def hist(data, wcs, imsz):
  coo = np.array(data, dtype='float64')[:,1:]
  foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
  H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
                            bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
  return H

if __name__ == '__main__':
  if False:
    skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    hdulist = pyfits.open('../data/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data

    step = 10
    offsets = np.load('../data/new/offsets300-1343_10_new.npy')
    length = offsets.shape[0]
    print length
    for i in range(length-1):
      offsets[i+1] = offsets[i+1] + offsets[i]
    offsets = np.append(offsets, np.array([[0.,0.]]), axis=0)
    print offsets.shape
    initial = int((co_data[300][0]-0.5)*1000)

    with open('../data/csv_list') as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          with open('../data/NUVphotons_5_cal_large_10_full.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile)
            data = []
            i = 1
            trange[0] = float(reader.next()[0])/1000
            id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
            for row in reader:
              time = int(row[0])
              if time/1000.>trange[0]+1342:
                print 'break'
                break
              row = cal_photon(offsets, initial, time, row, step)
              data.append(row[id])
              writer.writerow(row)
              if i%200000 == 0:
                H = hist(data, wcs, imsz)
                count += H
                trange[1] = float(data[-1][0])/1000
                data=H = None
                data = []
                gc.collect()
                print i
              i+=1
            if len(data)>0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_5_cal_large_10_sep_full.fits', clobber=False)


  if True:
    skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    skyrange = [24,15]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data

    step = 10
    offsets = np.load('../data/05_r/cata/offsets300-1343_10_new.npy')
    print offsets.shape
    initial = int((co_data[300][0]-0.5)*1000)

    with open('../data/csv_list05_r') as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          with open('../data/05_r/cata/NUVphotons_05_cal_cata_10_full.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile)
            data = []
            i = 1
            trange[0] = float(reader.next()[0])/1000
            id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
            for row in reader:
              time = int(row[0])
              if time/1000.>trange[0]+1342:
                print 'break'
                break
              row = cal_photon_r(offsets, initial, time, row, step)
              data.append(row[id])
              writer.writerow(row)
              if i%200000 == 0:
                H = hist(data, wcs, imsz)
                count += H
                trange[1] = float(data[-1][0])/1000
                data=H = None
                data = []
                gc.collect()
                print i
              i+=1
            if len(data)>0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_gPr_cata_10_sep_full_flip.fits', clobber=False)


  if False:
    skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    skyrange = [24,15]
    skypos, skyrange = ([0, 0], [1.33333336,1.33333336])
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.0016666667)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.0016666667)
    print imsz

    hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data

    initial = int((co_data[300][0]-0.5)*1000)

    with open('../data/csv_list05_device') as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file) 
          data = []
          i = 1
          trange[0] = float(reader.next()[0])/1000
          id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
          for row in reader:
            time = int(row[0])
            if time/1000.>trange[0]+1342:#780:
              print 'break'
              break
            data.append(np.array(row,dtype='float64')[id])
            if i%200000 == 0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
              data=H = None
              data = []
              gc.collect()
              print i
            i+=1
          if len(data)>0:
            H = hist(data, wcs, imsz)
            count += H
            trange[1] = float(data[-1][0])/1000
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_86_gPr_sep_device_r0_d_new.fits', clobber=False)
    #hdulist.writeto('../fits/count_map_05_gPr_sep_big.fits', clobber=False)

#focal plane 
  if False:
    skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    skyrange = [24,15]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
    co_data = hdulist[1].data

    step = 10
    offsets = np.load('../data/05_r/cata/offsets300-1343_10_new.npy')
    print offsets.shape
    initial = int((co_data[300][0]-0.5)*1000)


    skypos_f, skyrange_f = ([0, 0], [1.33333336,1.33333336])
    imsz_f = imagetools.deg2pix(skypos_f, skyrange_f, 0.0016666667)
    count_f = np.zeros(imsz_f)
    wcs_f = imagetools.define_wcs(skypos_f,skyrange_f,width=False,height=False,verbose=0,pixsz=0.0016666667)
    print imsz_f

    cata = spi.load_obj('stars05')
    cata_a = np.array(list(cata))
    cata_len = len(cata)
    cata_a = cata_a*np.pi/180.
    X, Y, Z = pycoo.spherical_to_cartesian(1, cata_a[:,1], cata_a[:,0])
    cata_car = np.array([X,Y,Z], dtype='float64').T
    rad = np.cos(0.00230*np.pi/180.)

    with open('../data/csv_list05_full') as f:
      csv_list = f.read().splitlines() 
      for csv_file in csv_list:
        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          data = []
          data_f = []
          i = 1
          trange[0] = float(reader.next()[0])/1000
          id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
          id_f = np.array([1,0,0,0,0,0,1,1,0,0,0])>0
          for row in reader:
            time = int(row[0])
            if time/1000.>trange[0]+1342:
              print 'break'
              break
            row = cal_photon_r(offsets, initial, time, row, step)
            if time<initial:
              pass
            else:
              center = pycoo.spherical_to_cartesian(1, row[9]*np.pi/180., row[8]*np.pi/180.)
              sep = np.dot(cata_car, center)
              coo = cata_a[sep>=rad,:]
              if len(coo) == 0:
                data_f.append(row[id_f])
            data.append(row[id])
            if i%200000 == 0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
              data=H = None
              data = []

              print 'select:%d'%len(data_f)
              if len(data_f) != 0:
                H_f = hist(data_f, wcs_f, imsz_f)
                count_f += H_f
                data_f=H_f=None
              data_f=[]
              gc.collect()
              print i
            i+=1
          if len(data)>0:
            H = hist(data, wcs, imsz)
            count += H
            trange[1] = float(data[-1][0])/1000

            H_f = hist(data_f, wcs_f, imsz_f)
            count_f += H_f
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_gPr_cata_10_sep_full_latest_new1342.fits', clobber=False)

    hdu = pyfits.PrimaryHDU(count_f)
    hdu = imagetools.fits_header('NUV', skypos_f, tranges, skyrange_f, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_gPr_full_cata_device_r0_latest_nostar_new1342.fits', clobber=False)


  if False:
    #skypos, skyrange = pos_range.get_pos_range(name_list='../data/name_list')
    skypos = [267.5, -29.0]
    skyrange = [26,24]
    tranges = []
    trange = [0, 0]
    imsz = imagetools.deg2pix(skypos, skyrange, 0.002)
    count = np.zeros(imsz)
    wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=0.002)
    print imsz

    step = 10
    scan_list = ['05','14','23','32','41','50','59','68']

    with open('../data/csv_list05_14_new') as f:
      csv_list = f.read().splitlines() 
      file_num = 0
      print file_num
      for csv_file in csv_list:
        scan_num = scan_list[file_num]
        print scan_num
        hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_000%s_0001-asprta.fits'%scan_num)
        offsets = np.load('../data/%s_new/cata/offsets1_10_new.npy'%scan_num)
        '''
        if file_num == 0:
          hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00005_0001-asprta.fits')
          offsets = np.load('../data/05_new/cata/offsets300-1343_10_new.npy')
          time_len = 1342
          print offsets.shape
        elif file_num == 1:
          hdulist = pyfits.open('../data/asprta/AIS_GAL_SCAN_00014_0001-asprta.fits')
          offsets = np.load('../data/14_new/cata/offsets300-1352_10_new.npy')
          time_len = 1351
          print offsets.shape
        '''
        co_data = hdulist[1].data
        initial = int((co_data[1][0]-0.5)*1000)

        print csv_file
        with open(csv_file, 'rb') as file:
          reader = csv.reader(file)
          data = []
          i = 1
          trange[0] = float(reader.next()[0])/1000
          id = np.array([1,0,0,0,0,0,0,0,1,1,0])>0
          for row in reader:
            time = int(row[0])
            '''
            if time/1000.>trange[0]+time_len:
              print 'break'
              break
            '''
            row = cal_photon_r(offsets, initial, time, row, step)
            data.append(row[id])
            if i%200000 == 0:
              H = hist(data, wcs, imsz)
              count += H
              trange[1] = float(data[-1][0])/1000
              data=H = None
              data = []
              gc.collect()
              print i
            i+=1
          if len(data)>0:
            H = hist(data, wcs, imsz)
            count += H
            trange[1] = float(data[-1][0])/1000
          file_num += 1
        tranges.append(trange)
        trange = [0, 0]

    print tranges
    print count.shape
    hdu = pyfits.PrimaryHDU(count)
    hdu = imagetools.fits_header('NUV', skypos, tranges, skyrange, hdu=hdu, wcs=wcs)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('../fits/count_map_05_68_gPr_cata_10_full_new_b.fits', clobber=False)

