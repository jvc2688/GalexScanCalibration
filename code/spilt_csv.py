import csv
import cPickle as pickle
import sys
import os

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

if __name__ == '__main__':
	csvname = sys.argv[1]

	output = "/beegfs/dw1519/galex/data/%s/split/lookup.pkl"%(csvname)
	dir = os.path.dirname(output)
	if not os.path.isfile(output):
		pass
	else:
		exit()
	if not os.path.exists(dir):
		os.makedirs(dir)
	else:
		print 'exists'
	look_up = {}
	i=0
	with open('/scratch/dw1519/galex/data/photon_list/%s.csv'%csvname, 'rb') as f:
		reader = csv.reader(f, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		filename = '%d.csv'%i
		w = open('/beegfs/dw1519/galex/data/%s/split/'%csvname+filename, 'wb')
		writer = csv.writer(w)
		row = reader.next()
		writer.writerow(row)
		last = int(row[0])
		look_up[row[0]] = i
		for row in reader:
			time = int(row[0])
			if time-last>=1000:
				i+=1
				filename = '%d.csv'%i
				last = time
				w.close()
				w = open('/beegfs/dw1519/galex/data/%s/split/'%csvname+filename, 'wb')
				writer = csv.writer(w)
				print time, filename
			writer.writerow(row)
			look_up[row[0]] = i
	save_obj(look_up, '/beegfs/dw1519/galex/data/%s/split/lookup'%csvname)
	os.remove('/scratch/dw1519/galex/data/photon_list/%s.csv'%csvname)






