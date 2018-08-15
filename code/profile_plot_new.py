import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import imagetools
from matplotlib.colors import LogNorm
import re
import glob
import os
import gc

def get_key(path):
	return int(re.split('_|/|\.', path)[-2])

def get_key(path):
	return int(re.split('/|\.', path)[-2])


def radial_profile1(data, center):
    r = np.sqrt((data[:,0])**2 + (data[:,1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r, minlength=4)
    radialprofile = tbin
    return radialprofile 

def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)
    print r
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

if __name__ == '__main__':
	if False:
		for i in range(100, 2900, 400):
			data_list=[]
			data_list1=[]
			for j in range(4):
				out = '/home/dw1519/dw1519/galex/plots/co239-10'
				data = np.load('{0}/{1}.npy'.format(out,i))/3600.
				data1 = np.load('{0}/ya/{1}.npy'.format(out,i))/3600.
				data_list.append(data)
				data_list1.append(data1)
				#info = np.load('/home/dw1519/dw1519/galex/plots/co32-10/info_{0}.npy'.format(i))
				#mask = info[:,1]>=2
			skypos = [0,0]
			skyrange = [0.02,0.02]
			pixsz = 0.0004
			imsz = imagetools.deg2pix(skypos, skyrange, pixsz).astype(int)
			wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(data[1:],1),1)
			H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
			                             bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
			#p = radial_profile(H, [99.5,99.5])
			p = radial_profile(H, [24.5,24.5])

			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(data1[1:],1),1)
			H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
			                             bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
			#p1 = radial_profile(H, [99.5,99.5])
			p1 = radial_profile(H, [24.5,24.5])
			x = np.arange(p.shape[0])*0.0004*3600
			p1 = p1/np.sum(p1[:15])
			p = p/np.sum(p[:15])
			plt.plot(x,p1, 'r', label='ya-correction')
			plt.plot(x,p, '--k', label='non-correction')
			plt.xlim(0,10)
			plt.xlabel('arcsec')
			plt.legend()
			plt.savefig('{0}/compare/p{1}.png'.format(out, i), dpi=190)
			plt.clf()

	if False:
		'''
		data_list = []
		data_list1 = []
		co_files =  sorted(glob.glob(out+"/[0-9]*.npy"), key=get_key)
		co_files1 =  sorted(glob.glob(out+"/cut/[0-9]*.npy"), key=get_key)
		print len(co_files)
		print len(co_files1)
		for cf, cf1 in zip(co_files, co_files1):
			co = np.load(cf)/3600.
			co1 = np.load(cf1)/3600.
			data_list.append(co)
			data_list1.append(co1)
		'''
		out = '/scratch/dw1519/galex/co/co239_1-10'
		path = os.path.dirname('{0}/compare/test'.format(out))
		if not os.path.exists(path):
			os.makedirs(path)
		co_files =  sorted(glob.glob(out+"/[0-9]*.npy"), key=get_key)
		num = get_key(co_files[-1])
		print num
		for i in range(100, 400, 400):
			data_list=[]
			data_list1=[]
			for j in range(4):
				data = np.load('{0}/cut/{1}.npy'.format(out,i+j*100))/3600.
				data1 = np.load('{0}/correct_new/{1}.npy'.format(out,i+j*100))/3600.
				data_list.append(data)
				data_list1.append(data1)
				#info = np.load('/home/dw1519/dw1519/galex/plots/co32-10/info_{0}.npy'.format(i))
				#mask = info[:,1]>=2]

			data = np.concatenate(data_list, axis=0)
			data1 = np.concatenate(data_list1, axis=0)
			a = 57./180*np.pi
			rot = np.array([[np.cos(a), -np.sin(a)],[np.sin(a), np.cos(a)]])
			#data = np.dot(data, rot)
			skypos = [0,0]
			skyrange = [0.04,0.04]
			pixsz = 0.0001
			imsz = imagetools.deg2pix(skypos, skyrange, pixsz).astype(int)
			wcs = imagetools.define_wcs(skypos,skyrange,width=False,height=False,verbose=0,pixsz=pixsz)
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(data[1:],1),1)
			H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
			                             bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
			p = radial_profile(H, [199.5,199.5])

			plt.imshow(H+10**-10,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[0.02*3600, -0.02*3600, 0.02*3600, -0.02*3600], origin='upper', norm=LogNorm(10**0,np.max(H)))  
			plt.colorbar()
			plt.xlim(-0.01*3600, 0.01*3600)
			plt.ylim(-0.01*3600, 0.01*3600)
			plt.ylabel(r'$\Delta gb$')
			plt.xlabel(r'$\Delta gl$')
			plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
			plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)
			#plt.plot(centroid[0]*3600,centroid[1]*3600,'+r', markersize=8)
			#plt.show()
			plt.title('t={0}-{1}s'.format(i/2-50, (i+300)/2))

			plt.savefig('{0}/compare/{1}.png'.format(out,i), dpi=190)
			plt.clf()


			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(data1[1:],1),1)
			H,xedges,yedges=np.histogram2d(foc[:,1]-0.5, foc[:,0]-0.5,\
			                             bins=imsz, range=([ [0,imsz[0]],[0,imsz[1]] ]))
			p1 = radial_profile(H, [199.5,199.5])

			plt.imshow(H+10**-10,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='equal', extent=[0.02*3600, -0.02*3600, 0.02*3600, -0.02*3600], origin='upper', norm=LogNorm(10**0,np.max(H)))  
			plt.colorbar()
			plt.xlim(-0.01*3600, 0.01*3600)
			plt.ylim(-0.01*3600, 0.01*3600)
			plt.ylabel(r'$\Delta gb$')
			plt.xlabel(r'$\Delta gl$')
			plt.axhline(y=0, color='b', linestyle='dashed', linewidth=1)
			plt.axvline(x=0, color='b', linestyle='dashed', linewidth=1)

			plt.title('t={0}-{1}s'.format(i/2-50, (i+300)/2))

			plt.savefig('{0}/compare/cut_{1}.png'.format(out,i), dpi=190)
			plt.clf()

			x = np.arange(p.shape[0])*0.0001*3600
			p1 = p1/np.sum(p1[:30])
			p = p/np.sum(p[:30])
			plt.plot(x,p1, 'r', label='correction')
			plt.plot(x,p, '--k', label='non-correction')
			plt.xlim(0,10)
			plt.xlabel('arcsec')
			plt.legend()
			plt.savefig('{0}/compare/p{1}.png'.format(out, i), dpi=190)
			plt.clf()

	if False:
		with open('../data/test_list') as f:
			name_list = f.read().splitlines()
		for name in name_list:
			scan = re.split('_', name)[3]
			qfile = '../co/qhist/'+scan+'.npy'
			if os.path.isfile(qfile):
				qhists = np.load(qfile)
			else:
				photon_list =  sorted(glob.glob("../data/"+name+"-cal-sec/split/*.csv"), key=get_key)
				#print photon_list
				qhists = []
				bins = np.arange(33)
				for f in photon_list:
					print f
					data = np.loadtxt(f,delimiter=',', usecols=[5])
					hist, edge = np.histogram(data, bins)
					qhists.append(hist)
				np.save('../co/qhist/'+scan+'.npy',qhists)
				qhists = np.array(qhists)
			print qhists.shape
			length=50
			num_seq = int(qhists.shape[0]/length)
			indices = (np.arange(num_seq)+1)*length
			s = np.split(qhists, indices, axis=0)
			s = np.array(map(lambda x: np.sum(x, axis=0), s)).T
			print s.shape
			nine=[]
			fif=[]
			ten=[]
			xl = s.T.shape[0]
			for i in range(s.T.shape[0]):
				sample=[]
				for j in range(32):
					sample += [j]*s[j,i]
				try:
					nine.append(np.percentile(sample,90))
					fif.append(np.percentile(sample,50))
					ten.append(np.percentile(sample,10))
				except IndexError:
					xl-=1
			s = s/(np.sum(s, axis=0).astype(float))
			plt.imshow(s,interpolation='None', cmap=plt.get_cmap('Greys'), aspect='auto',\
						 extent=[0, qhists.shape[0], 0, 32], origin='lower', norm=LogNorm())

			plt.plot((np.arange(xl, dtype=float)+0.5)*qhists.shape[0]/xl, nine,'r')
			plt.plot((np.arange(xl, dtype=float)+0.5)*qhists.shape[0]/xl, fif, 'r')
			plt.plot((np.arange(xl, dtype=float)+0.5)*qhists.shape[0]/xl, ten, 'r')
			plt.xlim(0,qhists.shape[0])
			plt.ylim(0,32)
			plt.title('scan {0}'.format(scan)) 
			plt.xlabel(r't [s]')
			plt.ylabel('q')
			#plt.savefig(out+'/xya_{0}.png'.format((i+1)*100), dpi=190)
			plt.savefig('../co/qhist/'+scan+'-norm_new.png', dpi=190)
			plt.clf()

	if True:
		with open('../co_pbs/np_new') as f:
			name_list = f.read().splitlines()
		for name in name_list:
			scan = re.split('_', name)[3]
			qfile = '../co/remove/'+scan+'.npy'
			if os.path.isfile(qfile):
				qhists = np.load(qfile)
			else:
				photon_list =  sorted(glob.glob("../data/"+name+"-cal-sec/split/*.csv"), key=get_key)
				#print photon_list
				qhists = []
				bins = np.arange(33)
				for f in photon_list:
					print f
					data = np.loadtxt(f,delimiter=',', usecols=[0,3,4,5])
					if len(data.shape)<2:
						continue
					mask = (data[:,2]>=2) & (data[:,3]>5)
					qhists.append([data[0,0], float(np.sum(mask))/data.shape[0]])
				np.save('../co/remove/'+scan+'.npy',qhists)
				qhists = np.array(qhists)
			print qhists.shape
			plt.plot((qhists[:,0]-qhists[0,0])/1000., qhists[:,1])
			plt.title('scan {0}'.format(scan)) 
			plt.xlabel(r't [s]')
			#plt.savefig(out+'/xya_{0}.png'.format((i+1)*100), dpi=190)
			plt.savefig('../co/remove/'+scan+'.pdf', dpi=190)
			plt.clf()
			qhists = None
			gc.collect()





