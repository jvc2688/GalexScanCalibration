import csv
import numpy as np

with open('../data/csv_list') as f:
  csv_list = f.read().splitlines() 
  for csv_file in csv_list:
    with open(csv_file, 'rb') as file:
      reader = csv.reader(file)
      initial_time = int(reader.next()[0])
      print initial_time
      for row in reader:
        time = int(row[0])
      print [initial_time, time]
      print (time - initial_time)