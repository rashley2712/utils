#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import trm.pgram

		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses matplotlib to plot light curves from rashley''s format CSV file.')
	parser.add_argument('datafile', nargs='+', type=str, help='Input data file(s)')
	
	 
	arg = parser.parse_args()
	print arg
	
	filename = arg.datafile[0]
	photometry = loadingSavingUtils.loadSingleChannelCSV(filename)
	
	""" Data is now loaded 
	"""

	# Fake some data
	x = numpy.linspace(0.,10.,100)
	f = 2.
	y = numpy.sin(2.*numpy.pi*f*x)
	e = 0.1*numpy.ones_like(x)
	
	ofac  = 300.
	hifac = 0.2
	
	x_values = []
	y_values = []
	y_errors = []
	for p in photometry:
		x_values.append(p['MJD'])
		y_values.append(p['magnitude'])
		y_errors.append(p['magnitudeError'])
		
	# Convert the x values into minutes
	x = numpy.array(x_values)
	minx = min(x) 
	x = (x - minx) * 1440
	
	y = numpy.array(y_values)
	e = numpy.array(y_errors)
	freq = numpy.arange(0, 1, 0.1)
	print freq
	f,p = trm.pgram.fwmls(x,y,e,ofac,hifac)
	#f,p = pgram.wmls(x,y,e,freq, ofac,hifac)
	
	print "Max red:",  max(p), " at freq: ", f[numpy.argmax(p)]
	
	matplotlib.pyplot.plot(f,p, 'r')
	
	
	
	matplotlib.pyplot.ylabel(r"P", size = 14)
	matplotlib.pyplot.xlabel(r"Frequency $(minutes^{-1})$", size = 14)
	
	
	filename = arg.datafile[1]
	photometry = loadingSavingUtils.loadSingleChannelCSV(filename)
	ofac  = 300.
	hifac = 0.1
	
	x_values = []
	y_values = []
	y_errors = []
	for p in photometry:
		x_values.append(p['MJD'])
		y_values.append(p['magnitude'])
		y_errors.append(p['magnitudeError'])
		
	# Convert the x values into minutes
	x = numpy.array(x_values)
	minx = min(x) 
	x = (x - minx) * 1440
	
	y = numpy.array(y_values)
	e = numpy.array(y_errors)
	freq = numpy.arange(0, 1, 0.1)
	print freq
	f,p = trm.pgram.fwmls(x,y,e,ofac,hifac)
	#f,p = pgram.wmls(x,y,e,freq, ofac,hifac)
	
	print "Max green:",  max(p), " at freq: ", f[numpy.argmax(p)]
	
	
	matplotlib.pyplot.plot(f,p, 'g')
	
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	
	fig.savefig('pgram.eps',dpi=100, format='eps')
	
	


