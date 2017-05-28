#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot
import scipy.signal as signal

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Attempt to calculate the probability of detecting RV variability based on the times of observations and a simulated period distribution. ')
	#parser.add_argument('filename', type=str, nargs='+', help='Filename of CSV file..')	
	parser.add_argument('--device', type=str, default='/xs', help = "Dump to the following device. Default is '/xs'.")
	parser.add_argument('--input', type=str, help='File containing a sample of periods to use as a period distribution.')
	
	arg = parser.parse_args()

	
	names = []
	samplePeriods = []
	if arg.input is not None:
		inputFile = open(arg.input, 'rt')
		for line in inputFile:
			values = line.strip().split('\t')
			names.append(values[0])
			samplePeriods.append(float(values[1]))
		inputFile.close()

	samplePeriods = [p/24. for p in samplePeriods]
	logSamplePeriods = numpy.log10(samplePeriods)
	print samplePeriods
	

	# the histogram of the data
	n, bins, patches = plt.hist(logSamplePeriods, 14, facecolor='green', alpha=0.75, cumulative=False)
	print n, bins, patches
	print numpy.sum(n), len(samplePeriods)
	plt.xlabel('$P_{orb}$ [d]')
	plt.ylabel('N')
	plt.title('Nebot period distribution')
	#plt.axis([0, 2, 0, 1])
	plt.grid(True)

	plt.show()


	periods = numpy.random.rand(1000)


	PGPlotWindow = ppgplot.pgopen(arg.device) 
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgslct(PGPlotWindow)   
	ppgplot.pgsci(1)
	ppgplot.pgask(False)
	ppgplot.pgsch(1.6)
	ppgplot.pgenv(0, 1, 0, 7, 0)
	ppgplot.pgsch(1.0)
	ppgplot.pghist(periods, 0, 1, 20, 2)
	ppgplot.pglab("P\dorb\u/days", "N", "")
	
	#ppgplot.pgclos()	
	
	
