import csv
import numpy
from itertools import islice
import pickle

def save_obj(obj, name ):
    with open('../data/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('../data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

if __name__ == '__main__':
	'''
	i = 0
	look_up = {}
	with open('../data/NUVphotons05.csv', 'rb') as f:
		reader = csv.reader(f)
		last = 0
		for row in reader:
			time = int(row[0])
			if time>last:
				look_up[row[0]] = i
				last = time
				print time,i
			i+=1
	save_obj(look_up, 'lookup')
	'''
	look_up = load_obj('lookup')
	print len(look_up.items())
	print look_up['1027851898990']
	print look_up['1027850555995']
