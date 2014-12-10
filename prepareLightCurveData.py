#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils

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
	
def binData(x, y, bin):
	newX = []
	newY = []
	
	newLength = len(x)/bin
	
	for i in range(newLength):
		xVal = 0
		yVal = 0
		for j in range(bin):
			xVal+= x[i*bin +j]
			yVal+= y[i*bin +j]
		xVal = xVal/bin
		yVal = yVal/bin
		newX.append(xVal)
		newY.append(yVal)
	
	return newX, newY
	
def filterData(dataArray):
	cleanData = []
	
	# This step removes the negative values
	for d in dataArray:
		items = len(d)
		if d[1]<0: continue
		if items>2:
			if d[2]<0: continue
		cleanData.append(d)
		
	return cleanData
	
def calcStats(data):
	np = numpy.array(data)
	values = np[:, 1]
	return (numpy.mean(values), numpy.std(values))
	
def appendPhotometry(existing, new):
	print "Appending photometry"
	colours = ['r', 'g', 'b']
	for c in colours:
		photometry = existing[c]
		newPhotometry = new[c]
		for n in newPhotometry:
			photometry.append(n)
			
	for c in colours:
		numDataPoints = len(existing[c])
		print c, "now has", numDataPoints, "data points"

	return existing
		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Reads CSV and FITS files and turns them into standard JSON files for loading into my various plotting programs,')
	parser.add_argument('datafile', nargs='+', type=str, help='Input data file(s)')
	parser.add_argument('-o', '--outfile', type=str, default='default', help='Name of the output file.')
	 
	arg = parser.parse_args()
	print arg
	colours = ['r', 'g', 'b']
	
	# Find last '.' in the filename and grab the file extension
	filename = str(arg.datafile[0])
	index = filename.rfind('.')
	extension = filename[index+1:].lower()
	filenameNoExtension = filename[:index]
	filetype = "unknown"
	if extension=="fits":
		filetype = 'fits'
	elif extension=="csv":
		filetype = "csv"
	
	if filetype == 'unknown':
		print "Sorry that filetype is unknown. Try CSV or FITS."
		
	
	if filetype == 'fits':
		print "Loading", filename
		photometry = loadingSavingUtils.loadFITSFile(filename)
		
		if (len(arg.datafile)>1):
			additionalFiles = arg.datafile[1:]
			for newFilename in additionalFiles:
				print "...also loading", newFilename
				additionalPhotometry = loadingSavingUtils.loadFITSFile(newFilename)
				
				photometry = appendPhotometry(photometry, additionalPhotometry)
	
	
	for c in colours:
		photometry[c], numRemoved = loadingSavingUtils.removeNegativeValues(photometry[c])
		if numRemoved>0: print "..removed %d negative values from %s."%(numRemoved, c)
	
	for c in colours:
		photometry[c], numRemoved = loadingSavingUtils.removeZeroValues(photometry[c])
		if numRemoved>0: print "..removed %d zero values from %s."%(numRemoved, c)
	
	
	if arg.outfile == 'default':
		outputFilename = filenameNoExtension + ".csv"
	else:
		outputFilename = arg.outfile
	loadingSavingUtils.writeCSV(outputFilename, photometry)
	
	
	
	sys.exit()
	
	
	for datafile in arg.datafiles:
	
		inputFile = astropy.io.fits.open(datafile)
		
		c = colours[0]
		
		headers = inputFile[CCDs[c]].header
		data = inputFile[CCDs[c]].data
		columns = inputFile[CCDs[c]].columns
		for item in data:
			time = item[columns.names.index("MJD")]
			reading = [time]			
			for col in fitsColumns:
				value = item[columns.names.index(col)]
				reading.append(value)
			if (arg.errors):
				for col in errorColumns:
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
			for col in fitsColumns:
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
			for col in fitsColumns:
				value = item[columns.names.index(col)]
				reading.append(value)
			try: 
				zeros = reading.index(0)
			except ValueError:
				blues.append(reading)
		
			
		inputFile.close()
	
	
	""" Data is now loaded 
	"""
	
		
	x_values = []
	y_values = []
	for r in reds:
		x_values.append(r[0])
		if len(fitsColumns)>1:
			if arg.m == True:
				y_values.append(-2.5 * math.log10(r[1]/r[2]))
			else:
				y_values.append(r[1]/r[2])
		else:
			y_values.append(r[1])
		if (arg.errors):
			errorValues = r[3]
		
	
	if (arg.bin!=1):
		x_values, y_values = binData(x_values, y_values, arg.bin)
	
	colourPlotData = {}
	colourPlotData['times'] = x_values
	colourPlotData['r'] = y_values
	
	
	MJDoffset = int(min(x_values))
	print "MJD offset:",MJDoffset
	
	if (arg.m == True):
		mean = numpy.mean(y_values)
		y_values = [y - mean for y in y_values]
	
	x_values = [x - MJDoffset for x in x_values]
	
	
	if (not arg.xb):
		matplotlib.pyplot.figure(figsize=(12, 12))
		axes = matplotlib.pyplot.subplot(3, 1, 3)
	else:
		matplotlib.pyplot.figure(figsize=(12, 8))
		axes = matplotlib.pyplot.subplot(2, 1, 2)
		
	#matplotlib.pyplot.plot(x_values, y_values, 'r.', label = 'i')
	matplotlib.pyplot.errorbar(x_values, y_values, yerr=0.3)
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
	ylabel_str = "$"
	if len(fitsColumns)==1: ylabel_str+= fitsColumns[0]
	else: ylabel_str+= fitsColumns[0] + " / " + fitsColumns[1]
	ylabel_str+= "$"
	matplotlib.pyplot.ylabel(r"Relative flux: $i$", size = 16)
	if (arg.m == True): 
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.ylabel(r"$i_{mag}$", size = 18)
	

	x_values = []
	y_values = []
	for g in greens:
		if len(fitsColumns)>1:
			if arg.m == True:
				if (g[1]<0): continue
				y_values.append(-2.5 * math.log10(g[1]/g[2]))
			else:
				y_values.append(g[1]/g[2])
		else:
			y_values.append(g[1])
		x_values.append(g[0])
		
	#y_values = astropy.stats.funcs.sigma_clip(y_values, sig=2, iters=1)
	if (arg.m == True):
		mean = numpy.mean(y_values)
		y_values = [y - mean for y in y_values]
	
	if (arg.bin!=1):
		x_values, y_values = binData(x_values, y_values, arg.bin)
	
	colourPlotData['g'] = y_values
	
		
	x_values = [x - MJDoffset for x in x_values]
	
	if (not arg.xb):
		axes = matplotlib.pyplot.subplot(3, 1, 2)
	else:
		axes = matplotlib.pyplot.subplot(2, 1, 1)
		
	matplotlib.pyplot.plot(x_values, y_values, 'g.', label = 'g')
	matplotlib.pyplot.ylabel(r"Relative flux: $g$", size = 16)
	if (arg.m == True): 
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.ylabel(r"$g_{mag}$", size = 18)
	
	blues = filterData(blues)
	
	if (not arg.xb):
		x_values = []
		y_values = []
		for b in blues:
			x_values.append(b[0])
			if len(fitsColumns)>1:
				if arg.m == True:
					y_values.append(-2.5 * math.log10(b[1]/b[2]))
				else:
					y_values.append(b[1]/b[2])
			else:
				y_values.append(b[1])
			
		#y_values = astropy.stats.funcs.sigma_clip(y_values, sig=2, iters=1)
		if (arg.m == True):
			mean = numpy.mean(y_values)
			y_values = [y - mean for y in y_values]
		
		if (arg.bin!=1):
			x_values, y_values = binData(x_values, y_values, arg.bin)
			
		colourPlotData['b'] = y_values
		colourPlotData['b_times'] = x_values
		
		x_values = [x - MJDoffset for x in x_values]
		
		matplotlib.pyplot.subplot(3, 1, 1)
		matplotlib.pyplot.plot(x_values, y_values, 'b.', label = 'u')
		matplotlib.pyplot.ylabel(r"Relative flux: $u$", size = 16)
		if (arg.m == True): 
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.ylabel(r"$u_{mag}$", size = 18)
		
	# Now check everything with the defaults:
	fig = matplotlib.pyplot.gcf()
	DPI = fig.get_dpi()
	print "DPI:", DPI
	DefaultSize = fig.get_size_inches()
	print "Default size in Inches", DefaultSize
	print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

	#matplotlib.pyplot.subplots_adjust(hspace=0.0, bottom=0.125)
	matplotlib.pyplot.show()
	
	fig.savefig('lightcurves.eps',dpi=100, format='eps')
	fig.savefig('lightcurves.png',dpi=100, format='png')
	
	# Now try a two-colour plot....
	x_values = colourPlotData['times']
	x_values = [x - MJDoffset for x in x_values]
	
	r_values = colourPlotData['r']
	g_values = colourPlotData['g']
	
	if (not arg.xb):
		b_values = colourPlotData['b']
		b_times = colourPlotData['b_times']
		
	col = zip(r_values, g_values)
	y_values = [n[1] - n[0] for n in col]
	
	if (arg.zero):
		mean = numpy.mean(y_values)
		y_values = [y-mean for y in y_values]
	
	if (not arg.xb):
		matplotlib.pyplot.figure(figsize=(12, 8))
		axes = matplotlib.pyplot.subplot(2, 1, 2)
	else:
		matplotlib.pyplot.figure(figsize=(12, 5))
		axes = matplotlib.pyplot.subplot(1, 1, 1)
	
	
	matplotlib.pyplot.plot(x_values, y_values, 'k.', label = 'g-i')
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.ylabel(r"$(g-i)_{mag}$", size = 18)
	
	
	if (not arg.xb):
		axes = matplotlib.pyplot.subplot(2, 1, 1, sharex = axes)
		
		col = zip(g_values, b_values)
		y_values = [n[1] - n[0] for n in col]
	
		if (arg.zero):
			mean = numpy.mean(y_values)
			y_values = [y-mean for y in y_values]
		
		x_values = [x - MJDoffset for x in b_times]
		
		matplotlib.pyplot.plot(x_values, y_values, 'k.', label = 'u-g')
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.ylabel(r"$(u-g)_{mag}$", size = 18)
		
	fig = matplotlib.pyplot.gcf()
	
	matplotlib.pyplot.show()
	fig.savefig('colourcurves.eps',dpi=100, format='eps')
	fig.savefig('colourcurves.png',dpi=100, format='png')
	
