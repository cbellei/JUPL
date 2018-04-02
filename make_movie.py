import csv
import matplotlib.pyplot as plt
import sys
import numpy as np
import os

# matplotlib.rc('font', size=24)
# matplotlib.rc('font', family='Arial')

try:
	nspec = int(sys.argv[1])
	maxn = sys.argv[2]
	dframe = sys.argv[3]
except:
	pass

prefix = "./output/"

def movie(nspec=4,maxn=10,dframe=1):
	"""nspec: number of ion species
	maxn: maximum number of files to scan through
	dframe: interval between file (for generating movie)"""

	param = list(csv.reader(open(prefix + 'parameters.csv')))[0]

	dtm = np.asarray(param[0]).astype(float)
	nz  = np.asarray(param[1]).astype(int)
	dt_print = np.asarray(param[2]).astype(float)
	tm_quiet = np.asarray(param[3]).astype(float)
	L = 1.e2 * np.asarray(param[4]).astype(float)  # in cm

	x = list( csv.reader(open(prefix + 'r.csv')) )
	x = np.asarray(x).astype(float)
	print "done with x"

	vel = dict()
	den = dict()
	temp = dict()

	for i in range(nspec+1):
		filename = prefix + "vel" + str(i+1) + ".csv"
		vel[i] = list( csv.reader(open(filename)) )
		vel[i] = np.asarray(vel[i]).astype(float)
		print "reading file " + filename

		filename = prefix + "den" + str(i+1) + ".csv"
		den[i] = list( csv.reader(open(filename)) )
		den[i] = np.asarray(den[i]).astype(float)
		print "reading file " + filename

		filename = prefix + "temp" + str(i+1) + ".csv"
		temp[i] = list( csv.reader(open(filename)) )
		temp[i] = np.asarray(temp[i]).astype(float)
		print "reading file " + filename


	efield = list(csv.reader(open(prefix + 'efield.csv')))
	efield = np.asarray(efield).astype(float)
	print "done with efield"

	xrange = [0.,L]
	Trange = [0., 2.]

	nfiles = np.int( x.size / nz )

	n_total = int(np.ceil(len(x)/float(nz)))

	print nz, n_total

	f1 = plt.figure(figsize=(10,10))
	f2 = plt.figure(figsize=(10,5))


	print "rmin = ", min(x)

	count = 0

	U = dict()
	N = dict()
	T = dict()
	lspec = ["C","D","3He","el"]
	clr = ["b","r","k"]
	if len(lspec)<nspec+1:
		print "---not enough labels - ABORTING"
		sys.exit(0)

	for n in range(1,min(int(maxn),n_total),int(dframe)):

		print n
		count += 1

		time = dt_print * (n-1)

		start = (n - 1) * nz; finish = start + nz
		mu = x[start:finish]

		for i in range(nspec+1):
			U[i] = vel[i][start:finish]
			N[i] = den[i][start:finish]
			T[i] = temp[i][start:finish]

		Efield = efield[start:finish]
#		p1 = n1 * 1.e6 * T1 * 1.6e-19
#		pe = ne * 1.e6 * Te * 1.6e-19
#  		Mu_ion = mu_ion[start:finish]

		mu = 1.e2 * mu

		plt.figure(1) #select figure 1

		plt.subplot(2,2,1)
		for i in range(nspec):
			plt.semilogy(mu,N[i]/1.e23,linewidth=2,label=lspec[i],color=clr[i])
		plt.xlabel('R [cm]')
		#ylabel('Density [1.e22 cm-3]')
		plt.ylabel('Density [1.e23 cm-3]')
		plt.xlim(xrange)
		plt.ylim([1.e-6, 5.e-3])
		title = 't = ' + str(1.e9*time) + ' ns'
		plt.legend(loc='lower right')
		plt.grid()
		plt.title(title)
		plt.locator_params(axis='x', nbins=4)


		plt.subplot(2,2,2)
		for i in range(nspec):
			plt.plot(mu,U[i]/1.e5,linewidth=2,color=clr[i])
		plt.xlabel('R [cm]')
		plt.ylabel('Velocity [10^7 cm/s]')
		plt.xlim(xrange)
		plt.ylim([-4,8])
		plt.grid()
		plt.locator_params(axis='x', nbins=4)

		plt.subplot(2,2,3)
		for i in range(nspec+1):
			plt.plot(mu,1.e-3*T[i],linewidth=2,color=clr[i])
		plt.xlabel('R [cm]')
		plt.ylabel('Temperature [keV]')
		plt.ylim(Trange)
		# ylim([5, 1.2*np.max(T2)])
		plt.xlim(xrange)
		plt.locator_params(axis='x', nbins=4)

		plt.subplot(2,2,4)
#		plt.plot(mu,Mu_ion/1.6e-19,linewidth=2)
		plt.plot(mu,Efield/1.e6,linewidth=2)
		plt.xlabel('R [cm]')
		plt.ylabel('E-field [MV/m]')
		plt.xlim(xrange)
		plt.locator_params(axis='x', nbins=4)

		filename = prefix + 'foo' + str(count).zfill(3) + '.jpg';

		print filename

		plt.savefig(filename)

		plt.clf()

	os.system("ffmpeg -y -i './output/foo%03d.jpg' output.m4v")
	#os.system("ffmpeg -y -i 'friction%03d.jpg' friction.m4v")
	#os.system("avconv -i 'foo%03d.jpg' -r 10 movie.avi")
	os.system("rm -f ./output/foo*.jpg")

	print "done"


if __name__ == "__main__":
	print "------------"
	print "nspec = " + str(nspec)
	print "maxn = " + maxn
	print "dframe = " + dframe
	print "------------"
	movie(nspec,maxn,dframe)
