#!/usr/bin/env python
import sys, os
import numpy, math, scipy.optimize
import argparse
import matplotlib.pyplot

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Compare the cumulative histograms of two period population samples.')
	parser.add_argument('file', type=str, help= 'File containing the period data (tsv).')
	args = parser.parse_args()
	print args

	periodFile = open(args.file, 'rt')

	names = []
	periods = []

	for index, f in enumerate(periodFile):
		params = f.strip().split('\t')
		try:
			name = str(params[0])
			period = float(params[1])
			names.append(name)
			periods.append(period)
		except:
			print "Invalid data on line", index
			continue

	print "Loaded %d periods."%(len(names))
	print

	periodFile.close()

	logPeriods = numpy.log10(periods)

	for name, period, logP in zip(names, periods, logPeriods):
		print(name, period, logP)


	# Draw the histogram of the input data
	figure1 = matplotlib.pyplot.figure(figsize=(7,10))
	binwidth=0.5
	data = logPeriods
	n, bins, patches = matplotlib.pyplot.hist(data, bins=numpy.arange(min(data), max(data) + binwidth, binwidth), facecolor='green', alpha=0.75, cumulative=False)
	print n, bins, patches
	print "Total in bins:", sum(n)
	matplotlib.pyplot.xlabel('$log_{10}(P_{orb})$ [d]',fontsize=18)
	matplotlib.pyplot.ylabel('N',fontsize=18)
	matplotlib.pyplot.title('Observed period distribution',fontsize=18)
	matplotlib.pyplot.grid(True)

	midpoints = []
	for index in range(len(bins)-1):
		midpoint = (bins[index] + bins[index+1]) / 2.0
		print index, bins[index], bins[index+1], midpoint
		midpoints.append(midpoint)

	matplotlib.pyplot.scatter(midpoints, n)

	def gaussian(x, a0, a1, a2):
		return a0  * numpy.exp(-.5 * ((x-a1)/a2)**2)

	a0 = max(n)
	a1 = numpy.mean(logPeriods)
	a2 = max(logPeriods) - min(logPeriods)

	print "Guess", a0, a1, a2


	guess = [a0, a1, a2]
	x = numpy.array(midpoints)
	y = numpy.array(n)
	try:
		results, covariance = scipy.optimize.curve_fit(gaussian, x, y, guess)
	except ValueError:
		print("Fit failed. Try to tweak")

	errors = numpy.sqrt(numpy.diag(covariance))
	(a0, a1, a2) = results

	xFit = numpy.arange(min(bins), max(bins), 0.1)
	yFit = gaussian(xFit, a0, a1, a2)

	matplotlib.pyplot.plot(xFit, yFit)

	print results

	print "logP median", numpy.median(logPeriods)
	print "logP mean", numpy.mean(logPeriods)

	matplotlib.pyplot.show()

	if args.save is not None:
		figure1.savefig(arg.save, bbox_inches='tight')
