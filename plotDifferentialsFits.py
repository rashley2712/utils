#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse
import astropy.io.fits

def findMatchingTime(data, target):
	reading = (0, -1)
	
	for d in data:
		time = d[0]
		if time==target:
			reading = (time, d[1])
			continue
	
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
	parser.add_argument('datafile', type=str, help='Input data file for object 1 (this time a FITS file.')
	arg = parser.parse_args()
	print arg
	
	fitsColumns = ['Counts_1', 'Counts_2']
	colours = ['r', 'g', 'b'];
	
	objects = [];
	
	for fits in fitsColumns:
		inputFile = astropy.io.fits.open(arg.datafile)
	
		objectPhotometry = {'r': None, 'g': None, 'b': None}
		
		headers = inputFile["CCD 1"].header
		data = inputFile["CCD 1"].data
		columns = inputFile["CCD 1"].columns
		reds = []
		for item in data:
			time = item[columns.names.index("MJD")]
			red = item[columns.names.index(fits)]
			if (red!=0):
				reds.append([time, red])
				
		headers = inputFile["CCD 2"].header
		data = inputFile["CCD 2"].data
		columns = inputFile["CCD 2"].columns
		greens = []
		for item in data:
			time = item[columns.names.index("MJD")]
			green = item[columns.names.index(fits)]
			if (green!=0):
				greens.append([time, green])
	
		headers = inputFile["CCD 3"].header
		data = inputFile["CCD 3"].data
		columns = inputFile["CCD 3"].columns
		blues = []
		for item in data:
			time = item[columns.names.index("MJD")]
			blue = item[columns.names.index(fits)]
			if (blue!=0):
				blues.append([time, blue])
	
		# Discard the first and last reading
		del reds[0]
		del reds[-1]
		del greens[0]
		del greens[-1]	
		del blues[0]
		del blues[-1]	
			
		objectPhotometry['r'] = reds
		objectPhotometry['g'] = greens
		objectPhotometry['b'] = blues
		objects.append(objectPhotometry)
			
		
		print fits + " - calcStats red", calcStats(reds)
		print fits + " - calcStats green", calcStats(greens)
		print fits + " - calcStats blue", calcStats(blues)
			

	inputFile.close()

	MJDoffset = int(reds[0][0])
	print "MJD offset:",MJDoffset
	
	
	x_values = []
	y_values = []
	times = []
	for r in objects[0]['r']:
		x_values.append(r[0] - MJDoffset)
		times.append(r[0])
		y_values.append(r[1])
		
	for index, mjd in enumerate(times):
		(closestTime, counts) = findMatchingTime(objects[1]['r'], mjd)
		if (counts == -1):
			print "No matching time"
		else:
			y_values[index] = y_values[index] / counts

	print "Diff mean [r]:" + str(numpy.mean(y_values)) + '[' + str(numpy.std(y_values)) + ']'

	matplotlib.pyplot.plot(x_values, y_values, 'r.', label = 'i')
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset))
	matplotlib.pyplot.ylabel('$Counts_{0} / Counts_{1}$')


	times = []
	x_values = []
	y_values = []
	for g in objects[0]['g']:
		x_values.append(g[0] - MJDoffset)
		times.append(g[0])
		y_values.append(g[1])
		
	for index, mjd in enumerate(times):
		(closestTime, counts) = findMatchingTime(objects[1]['g'], mjd)
		if (counts == -1):
			print "No matching time"
		else:
			y_values[index] = y_values[index] / counts

	print "Diff mean [g]:" + str(numpy.mean(y_values)) + '[' + str(numpy.std(y_values)) + ']'

	matplotlib.pyplot.plot(x_values, y_values, 'g.', label = 'g')
	
	times = []
	x_values = []
	y_values = []
	for b in objects[0]['b']:
		x_values.append(b[0] - MJDoffset)
		times.append(b[0])
		y_values.append(b[1])
		
	for index, mjd in enumerate(times):
		(closestTime, counts) = findMatchingTime(objects[1]['b'], mjd)
		if (counts == -1):
			print "No matching time"
		else:
			y_values[index] = y_values[index] / counts

	
	print "Diff mean [b]:" + str(numpy.mean(y_values)) + '[' + str(numpy.std(y_values)) + ']'
	matplotlib.pyplot.plot(x_values, y_values, 'b.', label = 'u')
	
	
	matplotlib.pyplot.legend()
	
	# Now check everything with the defaults:
	fig = matplotlib.pyplot.gcf()
	DPI = fig.get_dpi()
	print "DPI:", DPI
	DefaultSize = fig.get_size_inches()
	print "Default size in Inches", DefaultSize
	print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

	matplotlib.pyplot.show()
	
	fig.savefig('test2.eps',dpi=100, format='eps')
	
	
	
