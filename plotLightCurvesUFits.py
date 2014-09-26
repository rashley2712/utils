#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse
import astropy.io.fits
import astropy.stats

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

	parser = argparse.ArgumentParser(description='Uses matplotlib to plot light curves from the UCam log (that has been converted to fits format with "ulog2fits").')
	parser.add_argument('datafiles', nargs='+', type=str, help='Input data file(s)')
	parser.add_argument('-c', nargs='+', default = ['Counts_1'], type=str, help='Columns to plot')
	parser.add_argument('--channels', nargs='+', default = ['r', 'g', 'b'], type=str, help='Channels to plot')
	parser.add_argument('-m', action = 'store_true', help='Use magnitude scale')
	
	arg = parser.parse_args()
	print arg
	
	if len(arg.c)==0:
		fitsColumns = ["Counts_1"]
	else:
		fitsColumns = arg.c
		
	colours = ['r', 'g', 'b']
	CCDs = { 'r': 'CCD 1', 'g': 'CCD 2', 'b': 'CCD 3'}
	offSet = 0.1
	reds = []
	greens = []
	blues = []
	
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
		
	#y_values = astropy.stats.funcs.sigma_clip(y_values, sig=2, iters=1)
	
	MJDoffset = int(min(x_values))
	print "MJD offset:",MJDoffset
	
	if (arg.m == True):
		mean = numpy.mean(y_values)
		y_values = [y - mean for y in y_values]
	
	x_values = [x - MJDoffset for x in x_values]
	
	
	matplotlib.pyplot.figure(figsize=(12, 12))
	axes = matplotlib.pyplot.subplot(3, 1, 3)
	matplotlib.pyplot.plot(x_values, y_values, 'r.', label = 'i')
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset))
	ylabel_str = "$"
	if len(fitsColumns)==1: ylabel_str+= fitsColumns[0]
	else: ylabel_str+= fitsColumns[0] + " / " + fitsColumns[1]
	ylabel_str+= "$"
	matplotlib.pyplot.ylabel("$i_{mag}$")
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()

	x_values = []
	y_values = []
	for g in greens:
		x_values.append(g[0])
		if len(fitsColumns)>1:
			if arg.m == True:
				y_values.append(-2.5 * math.log10(g[1]/g[2]))
			else:
				y_values.append(g[1]/g[2])
		else:
			y_values.append(g[1])
		
	#y_values = astropy.stats.funcs.sigma_clip(y_values, sig=2, iters=1)
	if (arg.m == True):
		mean = numpy.mean(y_values)
		y_values = [y - mean for y in y_values]
		
	x_values = [x - MJDoffset for x in x_values]
	
	matplotlib.pyplot.subplot(3, 1, 2, sharex = axes)
	matplotlib.pyplot.plot(x_values, y_values, 'g.', label = 'g')
	matplotlib.pyplot.ylabel("$g_{mag}$")
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()

	
	blues = filterData(blues)
	
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
	
	x_values = [x - MJDoffset for x in x_values]
	
	matplotlib.pyplot.subplot(3, 1, 1)
	matplotlib.pyplot.plot(x_values, y_values, 'b.', label = 'u')
	matplotlib.pyplot.ylabel("$u_{mag}$")
	if (arg.m == True): matplotlib.pyplot.gca().invert_yaxis()



	
	# Now check everything with the defaults:
	fig = matplotlib.pyplot.gcf()
	DPI = fig.get_dpi()
	print "DPI:", DPI
	DefaultSize = fig.get_size_inches()
	print "Default size in Inches", DefaultSize
	print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

	#matplotlib.pyplot.subplots_adjust(hspace=0.0, bottom=0.125)
	matplotlib.pyplot.show()
	
	fig.savefig('test2.eps',dpi=100, format='eps')
	
	
	
