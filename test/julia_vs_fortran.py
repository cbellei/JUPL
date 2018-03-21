import numpy as np
import sys
import math
import csv


switch = 2

prefix1 = "./mION/"
prefix2 = "../output/"

if switch == 0:
	files = ["den1.csv", "den2.csv", "den3.csv",
				 "efield.csv", "r.csv",
				 "temp1.csv", "temp2.csv", "temp3.csv",
				 "vel1.csv", "vel2.csv", "vel3.csv"]
elif switch == 1:
	files = ["T.csv", "G1D.csv", "F1D.csv", "C1D.csv",
			 "U1D.csv", "U1D_p.csv", "U1D_c.csv"]
else:
	files = ["ne_cc.csv", "ni_cc.csv", "L_ab.csv", "L_ie.csv", "xiab.csv",
			 "taue.csv", "ke.csv", "nu_DT.csv", "k_DT.csv"]


def check_data(data1, data2):
	for i, (d1, d2) in enumerate(zip(data1, data2)):
		if not math.isclose(d1, d2, rel_tol=1e-5):
			print(suffix, i, d1, d2)
			sys.exit(0)

for suffix in files:
	print("suffix = " + suffix)
	file1 = prefix1 + suffix
	file2 = prefix2 + suffix

	with open(file1, 'r') as f1, open(file2, 'r') as f2:
		if switch==0:
			data1 = f1.read().splitlines()
			data1 = np.asarray(data1, float)
			data2 = f2.read().splitlines()
			data2 = np.asarray(data2, float)
			check_data(data1, data2)
		else:
			csv1 = csv.reader(f1)
			csv2 = csv.reader(f2)
			for i, (row1, row2) in enumerate(zip(csv1, csv2)):
				row1 = np.asarray(row1, dtype=float)
				row2 = np.asarray(row2, dtype=float)
				result = [math.isclose(r1, r2, rel_tol=1.e-4) for r1, r2 in zip(row1, row2)]
				if False in result and i>0:
					print(suffix, i+1)
					print(result)
					print("i = ", i+1)
					print("fortran")
					print(row1)
					print("julia")
					print(row2)
					sys.exit(0)
