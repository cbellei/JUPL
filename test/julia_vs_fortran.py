import numpy as np
import sys
import math
import csv


def compare(switch):

	prefix1 = "./mION/"
	prefix2 = "../output/"

	if switch == 0:
		files = ["den1.csv", "den2.csv", "den3.csv",
					 "r.csv", "temp1.csv", "temp2.csv", "temp3.csv",
					 "vel1.csv", "vel2.csv", "vel3.csv", "efield.csv"]
	elif switch == 1:
		files = ["F1D.csv", "C1D.csv",
				 "U1D.csv", "U1D_p.csv", "U1D_c.csv",
				 "T.csv", "G1D.csv"]
	else:
		files = ["ne_cc.csv", "ni_cc.csv", "Te_eV.csv", "T_eV.csv", "L_ab.csv", "L_ie.csv", "xiab.csv",
				 "taue.csv", "ke.csv", "nu_DT.csv", "k_DT.csv"]


	def check_data(data1, data2):
		for i, (row1, row2) in enumerate(zip(data1, data2)):
			row1 = np.asarray(row1, dtype=float)
			row2 = np.asarray(row2, dtype=float)
			result = [math.isclose(r1, r2, rel_tol=1.e-13) for r1, r2 in zip(row1, row2)]
			if False in result and i < 900:
				print(suffix, i+1)
				print(result)
				print("i = ", i+1)
				print("fortran")
				print(row1)
				print("julia")
				print(row2)
				sys.exit(0)

	for suffix in files:
		print("suffix = " + suffix)
		file1 = prefix1 + suffix
		file2 = prefix2 + suffix

		with open(file1, 'r') as f1, open(file2, 'r') as f2:
			csv1 = csv.reader(f1)
			csv2 = csv.reader(f2)
			check_data(csv1, csv2)


if __name__ == '__main__':
	switch = int(sys.argv[1])
	compare(switch)
