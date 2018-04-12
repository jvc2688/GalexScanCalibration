import os
import re
from os import listdir
from os.path import isfile, join

b1 = open('get_scan_raw6_batch_2', 'r')
name1 = set(b1.readlines()[0:-1])
b2 = open('get_scan_raw6_batch', 'r')
name2 = set(b2.readlines()[0:-1])

name1.update(name2)
name1 = list(name1)
print len(name1)

cwd = os.getcwd()
files = [f for f in listdir(cwd) if isfile(f)]
print len(files)
i=0
for name in name1:
  strings = re.split(' |/|\n', name)
  filename = strings[-2]
  if filename not in files:
    print '%s no'%filename
    os.system('scp broiler:/data2/galex/AIS_GAL_SCAN/raw6/%s .'%filename)
    i += 1
print i
