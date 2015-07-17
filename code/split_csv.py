import csv
import cPickle as pickle
import sys

def save_obj(obj, name ):
    with open('../data/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('../data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

if __name__ == '__main__':
	scan_num = sys.argv[1]
	look_up = {}
	i=0
	with open('../data/photon_list/NUVPhoton%s_full_new.csv'%scan_num, 'rb') as f:
		reader = csv.reader(f)
		filename = '%d.csv'%i
		w = open('../data/split%s_new/'%scan_num+filename, 'wb')
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
				w = open('../data/split%s_new/'%scan_num+filename, 'wb')
				writer = csv.writer(w)
				print time, filename
			writer.writerow(row)
			look_up[row[0]] = i
	save_obj(look_up, 'lookup%s_new'%scan_num)






