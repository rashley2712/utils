#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import trm.pgram as pgram

def findMatchingTime(data, target, tolerance):
	reading = (0, -1)
	distance = 1000
	for d in data:
		time = d[0]
		difference = abs(time - target)
		if difference<distance:
			distance = difference
			reading = (time, d[1])
			
	#print difference, tolerance
	#print "Target: ", target, "Found time: ", reading[0]
	if distance>tolerance:
		print "Target: ", target, "Found time: ", reading[0]
		reading = (target, -1)
		print "Panic!"
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
	parser.add_argument('datafile1', type=str, help='Input data file from new pipeline (a CSV file.')
	parser.add_argument('datafile2', type=str, help='Input data file from old pipeline (a FITS file.')
	parser.add_argument('-m', action = 'store_true', help='Use magnitude scale')
	
	arg = parser.parse_args()
	print arg


	colours = ['r', 'g', 'b']
	CCDs = { 'r': 'CCD 1', 'g': 'CCD 2', 'b': 'CCD 3'}
	
	objects = []
	
	file = arg.datafile1
	inputFile = open(file, 'r')
		
	objectPhotometry = {'r': None, 'g': None, 'b': None}
		
	reds = [];
	greens = [];
	blues = [];
		
	for n, line in enumerate(inputFile):
		values = line.split(',')
		if n==0:
			headings = values
		else:
			
			time = float(values[0])			
			try:
				red = float(values[1])
				reds.append([time, red])
			except:
				pass
	
			try:
				green = float(values[2])
				greens.append([time, green])
			except:
				pass
					
			try:
				blue = float(values[3])
				blues.append([time, blue])
			except:
				pass
	
	
	# Trim off the first and last data points
	reds.pop(0)
	reds.pop(-1)
	greens.pop(0)
	greens.pop(-1)
	blues.pop(0)
	blues.pop(-1)
	
	objectPhotometry['r'] = reds
	objectPhotometry['g'] = greens
	objectPhotometry['b'] = blues
	objects.append(objectPhotometry)
			
	inputFile.close()
	
	print str(file) + " - calcStats red", calcStats(reds)
	print str(file) + " - calcStats green", calcStats(greens)
	print str(file) + " - calcStats blue", calcStats(blues)
	
	
	objectPhotometry = {'r': None, 'g': None, 'b': None}
	reds = []
	greens = []
	blues = []	
	inputFile = astropy.io.fits.open(arg.datafile2)
		
	c = colours[0]
		
	headers = inputFile[CCDs[c]].header
	data = inputFile[CCDs[c]].data
	columns = inputFile[CCDs[c]].columns
	for item in data:
		time = item[columns.names.index("MJD")]
		reading = [time]			
		col = 'Counts_1'
		value = item[columns.names.index(col)]
		reading.append(value)
		try: 
			zeros = reading.index(0)
		except ValueError:
			reds.append(reading)
		
	c = colours[1]
	headers = inputFile[CCDs[c]].header
	data = inputFile[CCDs[c]].data
	columns = inputFile[CCDs[c]].columns
		
	for item in data:
		time = item[columns.names.index("MJD")]
		reading = [time]			
		col = "Counts_1"
		value = item[columns.names.index(col)]
		reading.append(value)
		try: 
			zeros = reading.index(0)
		except ValueError:
			greens.append(reading)
		
	c = colours[2]
	headers = inputFile[CCDs[c]].header
	data = inputFile[CCDs[c]].data
	columns = inputFile[CCDs[c]].columns
	for item in data:
		time = item[columns.names.index("MJD")]
		reading = [time]			
		col = "Counts_1"
		value = item[columns.names.index(col)]
		reading.append(value)
		try: 
			zeros = reading.index(0)
		except ValueError:
			blues.append(reading)
		

	objectPhotometry['r'] = reds
	objectPhotometry['g'] = greens
	objectPhotometry['b'] = blues
	objects.append(objectPhotometry)
	inputFile.close()
	
	print objects[0]['b']
	print objects[1]['b']
			
	MJDoffset = int(reds[0][0])
	print "MJD offset:",MJDoffset
	
	tolerance = abs(reds[10][0] - reds[9][0])
	
	matplotlib.pyplot.figure(figsize=(12, 12))
	axes = matplotlib.pyplot.subplot(3, 1, 3)
	
	
	x_values = []
	y_values = []
	times = []
	for r in objects[0]['r']:
		x_values.append(r[0] - MJDoffset)
		times.append(r[0])
		y_values.append(r[1])
		
	for index, mjd in enumerate(times):
		(closestTime, counts) = findMatchingTime(objects[1]['r'], mjd, tolerance)
		if (counts == -1):
			print "No matching time", mjd
		else:
			if (arg.m==True):
				y_values[index] = -2.5 * math.log10(y_values[index] / counts)
			else:
				y_values[index] = y_values[index] / counts
				
	if (arg.m == True):
		mean = numpy.mean(y_values)
		y_values = [y - mean for y in y_values]
	
	
	# Write to 'CSV'
	outputFile = open('out_r.csv', 'w')
	for index, x in enumerate(x_values):
		outputFile.write(str(x) + ", " + str(y_values[index]) + '\n')
	outputFile.close()

	print "Diff mean [r]:" + str(numpy.mean(y_values)) + '[' + str(numpy.std(y_values)) + ']'

	x_values_red = x_values
	y_values_red = y_values

	
	matplotlib.pyplot.plot(x_values, y_values, 'r.', label = 'i')
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.ylabel(r"Relative flux: $i$", size = 14)


	times = []
	x_values = []
	y_values = []
	for g in objects[0]['g']:
		x_values.append(g[0] - MJDoffset)
		times.append(g[0])
		y_values.append(g[1])
		
	for index, mjd in enumerate(times):
		(closestTime, counts) = findMatchingTime(objects[1]['g'], mjd, tolerance)
		if (counts == -1):
			print "No matching time"
		else:
			if (arg.m==True):
				y_values[index] = -2.5 * math.log10(y_values[index] / counts)
			else:
				y_values[index] = y_values[index] / counts

	if (arg.m == True):
		mean = numpy.mean(y_values)
		y_values = [y - mean for y in y_values]
	
	print "Diff mean [g]:" + str(numpy.mean(y_values)) + '[' + str(numpy.std(y_values)) + ']'

	x_values_green = x_values
	y_values_green = y_values

	matplotlib.pyplot.subplot(3, 1, 2, sharex = axes)
	matplotlib.pyplot.plot(x_values, y_values, 'g.', label = 'g')
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.ylabel(r"Relative flux: $g$", size = 14)
	
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
	
	
