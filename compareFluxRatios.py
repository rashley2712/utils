#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import loadingSavingUtils
import trm.pgram as pgram

def findMatchingTime(data, target):
	reading = (0, -1)
	distance = 1000
	for index, d in enumerate(data):
		time = d
		difference = abs(time - target)
		if difference<distance:
			distance = difference
			reading = (index, time)
			
	return reading

def findClosestTime(data, target):
	distance = 1000.
	reading = (0, 0)
	for d in data:
		time = d[0]
		gap = time - target
		if abs(gap) < distance:
			distance = abs(gap)
			reading = (time, d[1])
	return reading
	
def calcStats(data):
	np = numpy.array(data)
	values = np[:, 1]
	return (numpy.mean(values), numpy.std(values))
		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses matplotlib to makeshow the residuals of 2 light curves.')
	parser.add_argument('datafiles', nargs='+', type=str, help='Input data file(s)')
	parser.add_argument('-m', action = 'store_true', help='Use magnitude scale')
	parser.add_argument('--mjd', type=int, default = 0, help='Use this value as the MJD offset')
	
	arg = parser.parse_args()
	print arg


	colours = ['r', 'g', 'b']
	CCDs = { 'r': 'CCD 1', 'g': 'CCD 2', 'b': 'CCD 3'}
	
	autoPhotometry = loadingSavingUtils.loadSingleChannelCSV(arg.datafiles[0])
	tradPhotometry = loadingSavingUtils.loadSingleChannelCSV(arg.datafiles[1])
	
	tradTimes = []
	tradRatios = []
	for t in tradPhotometry:
		tradTimes.append(t['MJD'])
		tradRatios.append(t['fluxRatio'])
		
	
	x_values = []
	y_values = []
	differenceValues = []
	othery_values = []
	ratioValues = []
	for a in autoPhotometry:
		x_values.append(a['MJD'])
		y_values.append(a['fluxRatio'])
		index, match = findMatchingTime(tradTimes, a['MJD'])
		
		print match, a['MJD'], a['fluxRatio'], tradRatios[index]
		othery_values.append(tradRatios[index])
		differenceValues.append( a['fluxRatio'] - tradRatios[index] )
		ratioValues.append( a['fluxRatio'] / tradRatios[index] )
	
			
	MJDoffset = arg.mjd
	print "MJD offset:",MJDoffset
	
	
	matplotlib.pyplot.figure(figsize=(12, 12))
	axes = matplotlib.pyplot.subplot(3, 1, 3)
	
	
	matplotlib.pyplot.plot(x_values, ratioValues, 'r.', label = 'r')
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.ylabel(r"Auto F / Trad F", size = 14)
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
	

	autoPhotometry = loadingSavingUtils.loadSingleChannelCSV(arg.datafiles[2])
	tradPhotometry = loadingSavingUtils.loadSingleChannelCSV(arg.datafiles[3])
	
	tradTimes = []
	tradRatios = []
	for t in tradPhotometry:
		tradTimes.append(t['MJD'])
		tradRatios.append(t['fluxRatio'])
		
	
	x_values = []
	y_values = []
	differenceValues = []
	othery_values = []
	ratioValues = []
	for a in autoPhotometry:
		x_values.append(a['MJD'])
		y_values.append(a['fluxRatio'])
		index, match = findMatchingTime(tradTimes, a['MJD'])
		
		print match, a['MJD'], a['fluxRatio'], tradRatios[index]
		othery_values.append(tradRatios[index])
		differenceValues.append( a['fluxRatio'] - tradRatios[index] )
		ratioValues.append( a['fluxRatio'] / tradRatios[index] )
	
			
	MJDoffset = arg.mjd
	print "MJD offset:",MJDoffset
	
	
	axes = matplotlib.pyplot.subplot(3, 1, 2)
	
	
	matplotlib.pyplot.plot(x_values, ratioValues, 'g.', label = 'g')
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.ylabel(r"Auto F / Trad F", size = 14)


	autoPhotometry = loadingSavingUtils.loadSingleChannelCSV(arg.datafiles[4])
	tradPhotometry = loadingSavingUtils.loadSingleChannelCSV(arg.datafiles[5])
	
	tradTimes = []
	tradRatios = []
	for t in tradPhotometry:
		tradTimes.append(t['MJD'])
		tradRatios.append(t['fluxRatio'])
		
	
	x_values = []
	y_values = []
	differenceValues = []
	othery_values = []
	ratioValues = []
	for a in autoPhotometry:
		x_values.append(a['MJD'])
		y_values.append(a['fluxRatio'])
		index, match = findMatchingTime(tradTimes, a['MJD'])
		
		print match, a['MJD'], a['fluxRatio'], tradRatios[index]
		othery_values.append(tradRatios[index])
		differenceValues.append( a['fluxRatio'] - tradRatios[index] )
		ratioValues.append( a['fluxRatio'] / tradRatios[index] )
	
			
	MJDoffset = arg.mjd
	print "MJD offset:",MJDoffset
	
	
	axes = matplotlib.pyplot.subplot(3, 1, 1)
	
	
	matplotlib.pyplot.plot(x_values, ratioValues, 'b.', label = 'b')
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.ylabel(r"Auto F / Trad F", size = 14)

	
	matplotlib.pyplot.show()





	sys.exit()
	
	times = []
	x_values = []
	y_values = []
	for b in objects[0]['b']:
		x_values.append(b[0] - MJDoffset)
		times.append(b[0])
		y_values.append(b[1])
		
	for index, mjd in enumerate(times):
		(closestTime, counts) = findMatchingTime(objects[1]['b'], mjd, 3*tolerance)
		if (counts == -1):
			print "No matching time", mjd
		else:
			y_values[index] = y_values[index] / counts

	
	print "Diff mean [b]:" + str(numpy.mean(y_values)) + '[' + str(numpy.std(y_values)) + ']'
	matplotlib.pyplot.subplot(3, 1, 1, sharex = axes)
	matplotlib.pyplot.plot(x_values, y_values, 'b.', label = 'u')
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.ylabel(r"Relative flux: $u$", size = 14)
	
	
	
	# Now check everything with the defaults:
	fig = matplotlib.pyplot.gcf()
	DPI = fig.get_dpi()
	print "DPI:", DPI
	DefaultSize = fig.get_size_inches()
	print "Default size in Inches", DefaultSize
	print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

	matplotlib.pyplot.show()
	
	fig.savefig('test2.eps',dpi=100, format='eps')
	
	fig2 = matplotlib.pyplot.gcf()
	
	
	# Fake some data
	x = numpy.linspace(0.,10.,100)
	f = 2.
	y = numpy.sin(2.*numpy.pi*f*x)
	e = 0.1*numpy.ones_like(x)
	
	ofac  = 30.
	hifac = 0.1
	
	x_values = x_values_red
	y_values = y_values_red
	
	# Convert the x values into minutes
	x = numpy.array(x_values)
	minx = min(x) 
	x = (x - minx) * 1440
	
	y = numpy.array(y_values)
	e = numpy.ones_like(x)
	freq = numpy.arange(0, 1, 0.1)
	print freq
	f,p = pgram.fwmls(x,y,e,ofac,hifac)
	#f,p = pgram.wmls(x,y,e,freq, ofac,hifac)
	
	print "Max red:",  max(p), " at freq: ", f[numpy.argmax(p)]
	
	matplotlib.pyplot.plot(f,p, 'r')
	
	
	x_values = x_values_green
	y_values = y_values_green
	
	# Convert the x values into minutes
	x = numpy.array(x_values)
	minx = min(x) 
	x = (x - minx) * 1440
	
	y = numpy.array(y_values)
	e = numpy.ones_like(x)
	freq = numpy.arange(0, 1, 0.1)
	print freq
	f,p = pgram.fwmls(x,y,e,ofac,hifac)
	#f,p = pgram.wmls(x,y,e,freq, ofac,hifac)
	
	matplotlib.pyplot.plot(f,p, 'g')
	
	print "Max green:",  max(p), " at freq: ", f[numpy.argmax(p)]
	
	matplotlib.pyplot.show()
	
	fig2.savefig('pgram.eps',dpi=100, format='eps')
	
	
