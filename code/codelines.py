def convert(flat):
    r=[]
    mask = np.isfinite(flat)
    for i in range(800):
        rs = flat[i,:][mask[i,:]]
        if len(rs)>0:
            r.append(np.mean(rs))
        else:
            r.append(0)
    return np.array(r)

def convert1(flat):
    r=[]
    mask = np.isfinite(flat)
    for i in range(800):
        rs = flat[:,i][mask[:,i]]
        if len(rs)>0:
            r.append(np.mean(rs))
        else:
            r.append(0)
    return np.array(r)

r2 = convert(np.load('iterate2/flat8.npy'))

#plt.plot(r1/np.max(r1)/norm, '.', label='14.5-17')
#plt.plot(r3/np.max(r3)/norm, '.', label='15.5-17')
norm = ri/np.max(ri)
rs = hist/n
#mask = (x>380) & (x<420)
#plt.plot(x,ratio/np.mean(ratio[mask]),'.', label='scan/AIS')
plt.plot(rs/np.mean(rs[380:420]),'.', label='scan/AIS')
plt.plot(a,ri/np.max(ri)/norm, '.', label='pipeline X fraction')
plt.plot(a,rb/np.max(rb)/norm, '.', label='bkg')
plt.plot(a,r2/np.max(r2)/norm, '.', label='15-17')
plt.plot(a,rc3/np.max(rc3)/norm, '.', label='16-18')
plt.legend()
plt.ylim(0.4,1.2)
plt.show()

from scipy.ndimage.interpolation import rotate
roll=350.
a = np.arange(800)
stars = np.load('iterate3/stars_in.npy')
x = (stars[:,0]-5.)/1.33*800+400
ratio = np.power(10.,-stars[:,1]/2.5)
f,b = np.histogram(x,bins=800,range=(0,800),weights=ratio)
n,nb = np.histogram(x,bins=800,range=(0,800))
rs = f/n
flat = np.load('cut_iterate6/flat6.npy')
flat = np.swapaxes(flat,0,1)
flat = rotate(flat,-23.4,reshape=False,order=1,prefilter=False)
rtn = convert1(rotate(flat,-roll,reshape=False,order=1,prefilter=False))
flat = np.load('cut_iterate3/flat4.npy')
flat = np.swapaxes(flat,0,1)
flat = rotate(flat,-23.4,reshape=False,order=1,prefilter=False)
rt = convert1(rotate(flat,-roll,reshape=False,order=1,prefilter=False))
flat = np.load('iterate3/flat_in.npy')
flat = np.swapaxes(flat,0,1)
flat = rotate(flat,-23.4,reshape=False,order=1,prefilter=False)
rti = convert1(rotate(flat,-roll,reshape=False,order=1,prefilter=False))
flat = np.load('back_iterate1/flat9.npy')
flat = np.swapaxes(flat,0,1)
flat = rotate(flat,-23.4,reshape=False,order=1,prefilter=False)
rtb = convert1(rotate(flat,-roll,reshape=False,order=1,prefilter=False))
flat = np.load('iterate2/flat7.npy')
flat = np.swapaxes(flat,0,1)
flat = rotate(flat,-23.4,reshape=False,order=1,prefilter=False)
rt2 = convert1(rotate(flat,-roll,reshape=False,order=1,prefilter=False))
normt = rti/np.max(rti)
plt.plot(rs/np.mean(rs[380:420]),'.', label='scan/AIS')
plt.plot(a,rti/np.max(rti)/normt, '.', label='pipeline X fraction')
plt.plot(a,rtb/np.max(rtb)/normt, '.', label='bkg')
plt.plot(a,rt2/np.max(rt2)/normt, '.', label='15-17')
plt.plot(a,rt/np.max(rt)/normt, '.', label='16-18')
plt.plot(a,rtn/np.max(rtn)/normt, '.', label='15.5-17.5 big aperture')
plt.legend()
plt.ylim(0.4,1.2)
plt.show()


stars_all = np.load('/scratch/dw1519/galex/data/star_photon/cut_iterate6/stars_in_all.npy')
mask = stars_all[:,-1]>1.5
plt.scatter(stars_all[mask,0], stars_all[mask,-1],c=stars_all[mask,1], cmap=plt.get_cmap('jet'), alpha=0.8, edgecolors='face', s=0.5)
plt.ylabel('FWHM')
plt.xlabel('gl')
plt.ylim(0,8)
cb=plt.colorbar()
cb.set_label('NUV')
plt.show()
