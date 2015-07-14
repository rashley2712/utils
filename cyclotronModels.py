#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import numpy, math
import matplotlib.pyplot
import argparse
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot

def replaceExpChar(value):
	if 'D' in value:
		value = value.replace('D', 'E')
	
	return value
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads the output of a cyclotron model and plots it. ')
	parser.add_argument('inputfile', type=str, nargs='+', help='CSV style file containing the spectrum')
	 
	arg = parser.parse_args()
	print arg
	
	print "Astropy version:", astropy.__version__
	
	filename = arg.inputfile[0]
	
	inputfile = open(filename, 'r')
	freqs = []
	i_0s = []
	i_1s = []
	i_s = []
	for line in inputfile:
		if line[0] == '#':
			continue
		data = line.split()
		replaceExpChar(data[0])
		freqs.append(float(replaceExpChar(data[0])))
		i_0s.append(float(replaceExpChar(data[1])))
		i_1s.append(float(replaceExpChar(data[2])))
		i_s.append(float(replaceExpChar(data[3])))
		
		
	lowerFreq = min(freqs)
	upperFreq = max(freqs)
	upper_i0 = max(i_0s)
	lower_i0 = max(i_0s)
	upper_i1 = max(i_1s)
	lower_i1 = min(i_1s)
	upper_i = max(i_s)
	lower_i = min(i_s)
	
	mainPGPlotWindow = ppgplot.pgopen('/xs')	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgenv(lowerFreq, upperFreq, lower_i, upper_i, 0, 0)
	ppgplot.pgsci(1)
	
	ppgplot.pgline(freqs, i_s)
	ppgplot.pglab("frequency", "i_", "")
		
			
