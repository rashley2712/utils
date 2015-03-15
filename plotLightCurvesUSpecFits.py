#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import timeClasses


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
		

if __name__ == "__main__":
	print "Astropy version:", astropy.__version__


	parser = argparse.ArgumentParser(description='Uses matplotlib to plot light curves from the ULTRASPEC log (that has been converted to fits format with "ulog2fits").')
	parser.add_argument('datafiles', nargs='+', type=str, help='Input data file(s)')
	parser.add_argument('-c', nargs='+', default = ['Counts_1'], type=str, help='Columns to plot')
	parser.add_argument('-m', action = 'store_true', help='Use magnitude scale')
	parser.add_argument('--bin', type=int, default = 1, help='Binning factor')
	parser.add_argument('--zero', action = 'store_true', help='Remove the mean value from the plots.... Centering around zero.')
	parser.add_argument('--errors', action = 'store_true', help='Load and plot the error bars.')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	 
	arg = parser.parse_args()
	print arg
	
	if len(arg.c)==0:
		fitsColumns = ["Counts_1"]
	else:
		fitsColumns = arg.c
	
	if (arg.errors):
		errorColumns = ["Sigma_1", "Sigma_2"]
		
	if arg.e!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.e)
		print ephemeris
	else:
		hasEphemeris = False
		
	offSet = 0.1
	reds = []
	
	CCD = "CCD 1"
	
	sys.exit()
	
	for datafile in arg.datafiles:
	
		inputFile = astropy.io.fits.open(datafile)
		
		headers = inputFile['CCD 1'].header
		data = inputFile[CCD].data
		columns = inputFile[CCD].columns
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
	
	
	MJDoffset = int(min(x_values))
	print "MJD offset:",MJDoffset
	
	if (arg.m == True):
		mean = numpy.mean(y_values)
		y_values = [y - mean for y in y_values]
	
	x_values = [x - MJDoffset for x in x_values]
	
	
	matplotlib.pyplot.figure(figsize=(12, 8))
		
	matplotlib.pyplot.plot(x_values, y_values, 'r.')
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
	ylabel_str = "$"
	if len(fitsColumns)==1: ylabel_str+= fitsColumns[0]
	else: ylabel_str+= fitsColumns[0] + " / " + fitsColumns[1]
	ylabel_str+= "$"
	matplotlib.pyplot.ylabel(ylabel_str, size = 16)
	if (arg.m == True): 
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.ylabel(r"$i_{mag}$", size = 18)
	

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
		
		print "length g", len(g_values), "length b", len(b_values)
		col = zip(g_values, b_values)
		print col
		print "Final length", len(col)
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
	
