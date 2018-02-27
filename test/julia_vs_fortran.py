import numpy as np

files = ["den1.csv", "den2.csv", "den3.csv", 
			 "efield.csv", "r.csv", 
			 "temp1.csv", "temp2.csv", "temp3.csv",
			 "vel1.csv", "vel2.csv", "vel3.csv"]

prefix1 = "./mION/"
prefix2 = "../output/"

for suffix in files:
	print("suffix = " + suffix)
	file1 = prefix1 + suffix
	file2 = prefix2 + suffix
	
	with open(file1, 'r') as f1, open(file2, 'r') as f2:
		data1 = f1.read().splitlines()
		data1 = np.asarray(data1, float)
		data2 = f2.read().splitlines()
		data2 = np.asarray(data2, float)
	
	for i, (d1, d2) in enumerate(zip(data1, data2)):
		if(d1 != d2):
			print(suffix, i, d1, d2)