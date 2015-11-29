import csv
import numpy as np
import sys
import glob
from scipy import stats
import cPickle as pickle
import matplotlib.pyplot as plt
import spilt_csv as spi
import res_star
import gnomonic as gn
from astropy.io import fits as pyfits

def getKey(item):
	return item[2]

def getMag(item):
	return item[3]

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def running_sum(data, size, num):
	output = np.zeros(num)
	for i in range(num):
		output[i] = np.sum(data[i*size:(i+1)*size])
	return output

def cut_tail(data):
	length = data.shape[0]
	n = 1
	for i in range(1,length):
		if data[-i] == 0:
			n += 1
		else:
			break
	return np.min([length,length-n+2])

def likelihood(star_flux, res_flux):
	var = np.var(star_flux)
	return var, np.sqrt(np.mean((star_flux-res_flux)**2))

def main_0():
	name = 'AIS_GAL_SCAN_00032_0001'
	slope_list = []
	intercept_list = []
	star_list = []
	fout = open("../data/%s/star/extra/list_new/star_detector.csv"%name, 'ab')
	outwriter = csv.writer(fout)
	for i in range(5989):
		star_num = '%d'%i
		csv_file = "../data/%s/star/extra/list_new/star_%s-0.csv"%(name, star_num)
		star_csv = np.genfromtxt(csv_file, delimiter=',')
		x = np.array(star_csv[:,6], dtype='float64')
		y = np.array(star_csv[:,7], dtype='float64')
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
		slope_list.append(slope)
		intercept_list.append(intercept)
		star_list.append((i, slope, intercept))
		outwriter.writerow(['%d'%i, '%f'%slope, '%f'%intercept])
	save_obj(star_list, "../data/%s/star/extra/list_new/star_detector"%name)
	fout.close()

def main():

	name = 'AIS_GAL_SCAN_00014_0001'
	cata = spi.load_obj('../data/%s_starset_extra_full'%name)
	cata_a = np.array(list(cata))
	cata_len = len(cata)
	cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)
	star_list = load_obj("../data/%s/star/extra/list_new/star_detector"%name)
	star_list = sorted(star_list, key=getKey)
	print star_list
	plot_list = []
	stars = []
	istar = star_list[0]
	stars.append(istar[0])
	for i in range(1,len(star_list)):
		if star_list[i][2] - istar[2]<=100:
			stars.append(star_list[i][0])
		else:
			plot_list.append(stars)
			istar = star_list[i]
			stars = []
			stars.append(istar[0])

	group_num = 0
	for plot in plot_list:
		#f, axes = plt.subplots(len(plot), 1, squeeze=False)
		plot_num = 0
		xlim = [0,60] 
		if len(plot) == 1:
			continue
		flux_list = []
		for i in plot:
	#for i in range(397):
			star_count = []
			time_list = []
			bkg_time_list = []
			bkg_count = []
			star_num = '%d'%i
			csv_file = "../data/%s/star/extra/list_new/star_%s-0.csv"%(name, star_num)
			bkg_csv_file = "../data/%s/star/extra/bkg/star_bkg_%s-0.csv"%(name, star_num)
			star_file = open(csv_file, 'rb')
			bkg_file = open(bkg_csv_file, 'rb')
			star_reader = csv.reader(star_file)
			bkg_reader = csv.reader(bkg_file)
			last = 0
			for row in star_reader:
				time = float(row[0])#time = int(float(row[0])/1000.)
				if len(time_list) > 0:
					if time - time_list[0] > 59000:
						break
				if time - last > 1000:
					star_count.append(1)
					last = time
					time_list.append(last)
				else:
					star_count[-1] += 1
			last = 0
			bkg_count = np.zeros(len(star_count))
			for row in bkg_reader:
				time = float(row[0])#time = int(float(row[0])/1000.)
				for j in range(len(star_count)-1):
					if time>=time_list[j] and time<time_list[j+1]:
						bkg_count[j] += 1
						break
				if time>=time_list[-1]:
					bkg_count[-1] += 1
				'''
				if time - last > 1000:
					bkg_count.append(1)
					last = time
					bkg_time_list.append(last)
				else:
					bkg_count[-1] += 1
				'''

			star_flux = np.array(star_count) - np.array(bkg_count)*64./40.
			star_flux[star_flux<0] = 0
			if np.mean(star_flux)<5:
				continue
			#without bkg sub
			star_flux = np.array(star_count)
			time_array = (np.array(time_list)-time_list[0])/1000.
			star_file.close()
			bkg_file.close()
			'''
			plt.plot(time_array, star_flux)
			plt.xlabel('t/s')
			plt.ylabel('Number of Photons')
			plt.savefig('../plots/%s/star/star_flux/%s-0.png'%(name, star_num), dpi=190)
			plt.clf()
			'''
			flux_list.append((time_array, star_flux, i, cata_list[i, 5]))
		if len(flux_list)>1:
			flux_list = sorted(flux_list, key=getMag)
			f, axes = plt.subplots(len(flux_list), 1, squeeze=False)
			for plot_num in range(len(flux_list)):
				star_flux = flux_list[plot_num][1]
				time_array = flux_list[plot_num][0]
				axes[plot_num,0].plot(time_array, star_flux)
				axes[plot_num,0].text(0.95, 0.85, 'NUV: %.3f'%flux_list[plot_num][3], verticalalignment='bottom', horizontalalignment='right', transform=axes[plot_num,0].transAxes, color='green', fontsize=10)
				#axes[plot_num,0].set_ylabel('Num of Photons')
				ylim = axes[plot_num,0].get_ylim()
				if plot_num == 0:
					xlim = axes[plot_num,0].get_xlim()
				axes[plot_num,0].set_xlim([0, xlim[1]])
				axes[plot_num,0].set_ylim([ylim[0]+np.absolute(ylim[1]-ylim[0])*0.01, ylim[1]])
				if plot_num<len(flux_list)-1:
					plt.setp( axes[plot_num,0].get_xticklabels(), visible=False)
				else:
					axes[plot_num,0].set_xlabel('time/s')
				#plot_num+=1
			plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
			            wspace=0, hspace=0)
			plt.legend()
			plt.savefig('../plots/%s/star/extra/cons_no/cons%d_track.png'%(name, group_num),dpi=190)  
			plt.clf()
		group_num += 1

def main_new():

	name = 'AIS_GAL_SCAN_00014_0001'
	cata = spi.load_obj('../data/%s_starset_extra_full'%name)
	cata_a = np.array(list(cata))
	cata_len = len(cata)
	cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)
	star_list = load_obj("../data/%s/star/extra/list_new/star_detector"%name)
	star_list = sorted(star_list, key=getKey)
	print star_list
	plot_list = []
	stars = []
	istar = star_list[0]
	stars.append(istar[0])
	for i in range(1,len(star_list)):
		if star_list[i][2] - istar[2]<=100:
			stars.append(star_list[i][0])
		else:
			plot_list.append(stars)
			istar = star_list[i]
			stars = []
			stars.append(istar[0])

	group_num = 0
	for plot in plot_list:
		#f, axes = plt.subplots(len(plot), 1, squeeze=False)
		plot_num = 0
		xlim = [0,60] 
		if len(plot) == 1:
			continue
		flux_list = []
		for i in plot:
	#for i in range(397):
			star_count = []
			time_list = []
			bkg_time_list = []
			bkg_count = []
			star_num = '%d'%i
			csv_file = "../data/%s/star/extra/list_new/star_%s-0.csv"%(name, star_num)
			bkg_csv_file = "../data/%s/star/extra/bkg/star_bkg_%s-0.csv"%(name, star_num)

			star_data = np.genfromtxt(csv_file, delimiter=',')
			bkg_data = np.genfromtxt(bkg_csv_file, delimiter=',')

			star_count, edges = np.histogram(star_data[:,0]-star_data[0,0], bins=np.arange(62)*1000)

			bkg_count, edges = np.histogram(bkg_data[:,0]-star_data[0,0], bins=np.arange(62)*1000)

			star_flux = np.array(star_count) - np.array(bkg_count)*64./40.
			star_flux[star_flux<0] = 0
			if np.mean(star_flux)<5:
				continue
			time_array = np.arange(61)
			'''
			plt.plot(time_array, star_flux)
			plt.xlabel('t/s')
			plt.ylabel('Number of Photons')
			plt.savefig('../plots/%s/star/star_flux/%s-0.png'%(name, star_num), dpi=190)
			plt.clf()
			'''
			flux_list.append((time_array, star_flux, i, cata_list[i, 5]))
		if len(flux_list)>1:
			flux_list = sorted(flux_list, key=getMag)
			f, axes = plt.subplots(len(flux_list), 1, squeeze=False)
			for plot_num in range(len(flux_list)):
				star_flux = flux_list[plot_num][1]
				time_array = flux_list[plot_num][0]
				axes[plot_num,0].plot(time_array, star_flux)
				axes[plot_num,0].text(0.95, 0.85, 'NUV: %.3f'%flux_list[plot_num][3], verticalalignment='bottom', horizontalalignment='right', transform=axes[plot_num,0].transAxes, color='green', fontsize=10)
				#axes[plot_num,0].set_ylabel('Num of Photons')
				ylim = axes[plot_num,0].get_ylim()
				if plot_num == 0:
					xlim = axes[plot_num,0].get_xlim()
				axes[plot_num,0].set_xlim([0, xlim[1]])
				axes[plot_num,0].set_ylim([ylim[0]+np.absolute(ylim[1]-ylim[0])*0.01, ylim[1]])
				if plot_num<len(flux_list)-1:
					plt.setp( axes[plot_num,0].get_xticklabels(), visible=False)
				else:
					axes[plot_num,0].set_xlabel('time/s')
				#plot_num+=1
			plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
			            wspace=0, hspace=0)
			plt.legend()
			plt.savefig('../plots/%s/star/extra/cons_new/cons%d_track.png'%(name, group_num),dpi=190)  
			plt.clf()
		group_num += 1

def main_res():

	name = 'AIS_GAL_SCAN_00032_0001'

	asp = np.load('../data/photon_list/AIS_GAL_SCAN_00032_0001_asp_cal.npy')
	asp_uni, uni_index=np.unique(asp[:,0], True)
	asp_uni = asp[uni_index]
	ix_tmp = (np.round(asp_uni[:,0]-asp_uni[0,0])+1).astype(int)
	dead_tmp = np.zeros(ix_tmp.shape[0])
	scst_tmp = pyfits.open('../AIS_GAL_SCAN/scst/%s-scst.fits'%name)[1].data
	limit = scst_tmp['t_dead_nuv'].shape[0]-1
	ix_tmp[ix_tmp>limit] = limit
	dead_tmp = scst_tmp['t_dead_nuv'][ix_tmp]

	cata = spi.load_obj('../data/%s_starset_extra_full'%name)
	cata_a = np.array(list(cata))
	cata_len = len(cata)
	cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)
	star_list = load_obj("../data/%s/star/extra/list_new/star_detector"%name)
	star_list = sorted(star_list, key=getKey)
	#print star_list
	plot_list = []
	stars = []
	istar = star_list[0]
	stars.append(istar[0])
	for i in range(1,len(star_list)):
		if star_list[i][2] - istar[2]<=100:
			stars.append(star_list[i][0])
		else:
			plot_list.append(stars)
			istar = star_list[i]
			stars = []
			stars.append(istar[0])

	group_num = 0
	flux_stat = []
	for plot in plot_list:
		#f, axes = plt.subplots(len(plot), 1, squeeze=False)
		print('group_num:%d'%group_num)
		plot_num = 0
		xlim = [0,60] 
		if len(plot) == 1:
			group_num += 1
			continue
		flux_list = []
		for i in plot:
	#for i in range(397):
			star_count = []
			time_list = []
			bkg_time_list = []
			bkg_count = []
			star_num = '%d'%i
			csv_file = "../data/%s/star/extra/list_new/star_%s-0.csv"%(name, star_num)
			bkg_csv_file = "../data/%s/star/extra/bkg/star_bkg_%s-0.csv"%(name, star_num)
			star_co = cata_a[i, 0:2]

			star_data = np.genfromtxt(csv_file, delimiter=',')
			bkg_data = np.genfromtxt(bkg_csv_file, delimiter=',')

			star_count, edges = np.histogram(star_data[:,0]-star_data[0,0], bins=np.arange(62)*1000)

			bkg_count, edges = np.histogram(bkg_data[:,0]-star_data[0,0], bins=np.arange(62)*1000)

			index_0 = (star_data[0,0]-asp_uni[0,0]*1000)/5
			index_1 = (star_data[-1,0]-asp_uni[0,0]*1000)/5

			#print i, cata_list[i, 0]
			#print star_co
			#print np.median(star_data[:,-3:-1], axis=0)

			if index_1-index_0>12400:
				continue

			star_track = asp_uni[index_0:index_1+1,:]
			dead_t = dead_tmp[index_0:index_1+1]
			ra = np.ones(star_track.shape[0])*star_co[0]
			dec = np.ones(star_track.shape[0])*star_co[1]
			xi, eta = gn.gnomfwd_simple(ra, dec, star_track[:,1], star_track[:,2], -star_track[:,3], 1/36000., 0.0)
			res_list = res_star.get_res_new(xi, eta, '../fits/flat/NUV_flat.fits')
			if res_list == None:
				continue
			res_list_d = res_list*(1-dead_t)

			res_flux = running_sum(res_list_d, 200, 61)
			res_len = cut_tail(res_flux)

			star_flux = np.array(star_count) - np.array(bkg_count)*64./40.
			star_flux[star_flux<0] = 0
			star_len = cut_tail(star_flux)

			delta = star_len-res_len
			if delta>0 and delta<3:
				star_flux = star_flux[:star_len]
				res_flux = res_flux[:star_len]
			elif delta >= 3 or delta<=-3:
				continue
			else:
				star_flux = star_flux[:res_len]
				res_flux = res_flux[:res_len]

			if np.mean(star_flux)<5:
				continue

			res_flux = res_flux/np.mean(res_flux)
			star_flux = star_flux/np.mean(star_flux)
			time_array = np.arange(res_flux.shape[0])

			var, likeli = likelihood(star_flux, res_flux)

			mag = cata_list[i, 5]
			flux_stat.append((mag, likeli, var, star_flux, res_flux))
			print '%.3f  %.3f  %.3f'%(mag, var, likeli)
			'''
			plt.plot(time_array, star_flux)
			plt.xlabel('t/s')
			plt.ylabel('Number of Photons')
			plt.savefig('../plots/%s/star/star_flux/%s-0.png'%(name, star_num), dpi=190)
			plt.clf()
			'''
			flux_list.append((time_array, star_flux, i, cata_list[i, 5], res_flux))
		if len(flux_list)>1:
			flux_list = sorted(flux_list, key=getMag)
			f, axes = plt.subplots(len(flux_list), 1, squeeze=False)
			for plot_num in range(len(flux_list)):
				star_flux = flux_list[plot_num][1]
				time_array = flux_list[plot_num][0]
				res_flux = flux_list[plot_num][-1]
				axes[plot_num,0].plot(time_array, star_flux, '-o')
				axes[plot_num,0].plot(time_array, res_flux, '-o')
				axes[plot_num,0].text(0.95, 0.8, 'NUV: %.3f'%flux_list[plot_num][3], verticalalignment='bottom', horizontalalignment='right', transform=axes[plot_num,0].transAxes, color='green', fontsize=10)
				#axes[plot_num,0].set_ylabel('Num of Photons')
				ylim = axes[plot_num,0].get_ylim()
				if plot_num == 0:
					xlim = axes[plot_num,0].get_xlim()
				axes[plot_num,0].set_xlim([0, 61])
				axes[plot_num,0].set_ylim([ylim[0]+np.absolute(ylim[1]-ylim[0])*0.01, ylim[1]])
				if plot_num<len(flux_list)-1:
					plt.setp( axes[plot_num,0].get_xticklabels(), visible=False)
				else:
					axes[plot_num,0].set_xlabel('time/s')
				plt.setp( axes[plot_num,0].get_yticklabels(), visible=False)
			plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
			            wspace=0, hspace=0)
			plt.legend()
			plt.savefig('../plots/%s/star/extra/cons_res_flat/cons%d_track.png'%(name, group_num),dpi=190)  
			plt.clf()
		group_num += 1
		save_obj(flux_stat, '../plots/%s/star/extra/cons_res_flat/flux_stat'%(name))


def main_rms():
	name = 'AIS_GAL_SCAN_00014_0001'
	stat = load_obj('../plots/%s/star/extra/cons_res_cut/flux_stat'%(name))
	mag = []
	rms = []
	star = []
	for i in range(0, len(stat)):
		mag.append(stat[i][0])
		rms.append(stat[i][1])

	plt.plot(mag, rms, '.k')
	plt.xlabel('NUV Magnitude')
	plt.ylabel('RMS Error')
	plt.savefig('../plots/%s/star/extra/cons_res_cut/rms.png'%name, dpi=190)

def main_single():

	name = 'AIS_GAL_SCAN_00014_0001'

	cata = spi.load_obj('../data/%s_starset_extra_full'%name)
	cata_a = np.array(list(cata))
	cata_len = len(cata)
	cata_list = np.insert(cata_a, 0, np.arange(cata_len), axis=1)
	star_list = load_obj("../data/%s/star/extra/list_new/star_detector"%name)
	star_list = sorted(star_list, key=getKey)
	#print star_list
	plot_list = []
	stars = []
	istar = star_list[0]
	stars.append(istar[0])

	for i in range(5989):
		flux_list = []
		star_count = []
		time_list = []
		bkg_time_list = []
		bkg_count = []
		star_num = '%d'%i
		csv_file = "../data/%s/star/extra/list_new/star_%s-0.csv"%(name, star_num)
		bkg_csv_file = "../data/%s/star/extra/bkg/star_bkg_%s-0.csv"%(name, star_num)
		star_co = cata_a[i, 0:2]

		star_data = np.genfromtxt(csv_file, delimiter=',')
		bkg_data = np.genfromtxt(bkg_csv_file, delimiter=',')

		star_count, edges = np.histogram(star_data[:,0]-star_data[0,0], bins=np.arange(62)*1000)

		bkg_count, edges = np.histogram(bkg_data[:,0]-star_data[0,0], bins=np.arange(62)*1000)

		star_flux = np.array(star_count) - np.array(bkg_count)*64./40.
		#star_flux[star_flux<0] = 0
		star_len = cut_tail(star_flux)

		time_array = np.arange(star_flux.shape[0])

		mag = cata_list[i, 5]

		flux_list.append((time_array, star_flux, i, cata_list[i, 5]))
		flux_list = sorted(flux_list, key=getMag)
		f, axes = plt.subplots(len(flux_list), 1, squeeze=False)
		for plot_num in range(len(flux_list)):
			star_flux = flux_list[plot_num][1]
			time_array = flux_list[plot_num][0]
			res_flux = flux_list[plot_num][-1]
			axes[plot_num,0].plot(time_array, star_flux, '-o')
			axes[plot_num,0].text(0.95, 0.85, 'NUV: %.3f'%flux_list[plot_num][3], verticalalignment='bottom', horizontalalignment='right', transform=axes[plot_num,0].transAxes, color='green', fontsize=10)
			#axes[plot_num,0].set_ylabel('Num of Photons')
			ylim = axes[plot_num,0].get_ylim()
			if plot_num == 0:
				xlim = axes[plot_num,0].get_xlim()
			axes[plot_num,0].set_xlim([0, 61])
			axes[plot_num,0].set_ylim([ylim[0]+np.absolute(ylim[1]-ylim[0])*0.01, ylim[1]])

			axes[plot_num,0].set_xlabel('time/s')
			plt.legend()
			plt.savefig('../plots/%s/star/extra/single/star%d.png'%(name, flux_list[plot_num][2]),dpi=190)  
			plt.clf()

if __name__ == '__main__':
	main_res()
