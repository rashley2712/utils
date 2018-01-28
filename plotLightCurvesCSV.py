#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import timeClasses



if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses matplotlib to plot light curves from rashley''s format CSV file.')
	parser.add_argument('inputfiles', type=str, nargs='+', help='Input data in CSV format')
	parser.add_argument('--bin', type=int, default = 1, help='Binning factor')
	parser.add_argument('--zero', action = 'store_true', help='Remove the mean value from the plots.... Centering around zero.')
	parser.add_argument('--errors', action = 'store_true', help='Load and plot the error bars.')
	parser.add_argument('--phaseplot', action = 'store_true', help = 'Do a phased plot')
	parser.add_argument('-e', '--ephemeris', type=str, help='Optional ephemeris file')
	parser.add_argument('-y', '--invert', action = "store_true", help = 'Invert the y-axis (for magnitude plots).')
	parser.add_argument('--yvalue', type=str, help = "Specify the name of the y-axis plot value. Taken from the first row in th .csv file.")
	parser.add_argument('--yerror', type=str, help = "Specify the name of the y-axis error value. Taken from the first row in th .csv file.")

	arg = parser.parse_args()
	print arg

	if arg.ephemeris!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.ephemeris)
		print ephemeris
	else:
		hasEphemeris = False

	allData = []

	numPlots = 0
	for filename in arg.inputfiles:
		columnNames, photometry = loadingSavingUtils.loadNewCSV(filename)
		allData.append(photometry)
		numPlots+= 1


	columnIndex = 1
	errorIndex = 2
	if arg.yvalue:
		for index, header in enumerate(columnNames):
			if arg.yvalue==header:
				columnIndex = index
	if arg.yerror:
		for index, header in enumerate(columnNames):
			if arg.yerror==header:
				errorIndex = index
	""" Data is now loaded
	"""
	c = 'k'
	sizePerPlot = 4

	xColumn = columnNames[0]
	yColumn = columnNames[columnIndex]
	yErrors = columnNames[errorIndex]

	print allData

	matplotlib.pyplot.figure(figsize=(12, 12/1.6))

	for index, photometry in enumerate(allData):
		subPlot = index+1
		x_temp = photometry[xColumn]
		y_temp = photometry[yColumn]
		e_temp = photometry[yErrors]
		# remove empty values (with a value of '-1')
		print [(x, y, e) for (x, y, e) in zip(x_temp, y_temp, e_temp) ]
		x_values = []
		y_values = []
		y_errors = []
		for (x, y, e) in zip(x_temp, y_temp, e_temp):
			if y!=-1:
				x_values.append(x)
				y_values.append(y)
				y_errors.append(e)
		axes = matplotlib.pyplot.subplot(numPlots, 1, subPlot)
		print "subplot", numPlots, 1, subPlot

		matplotlib.pyplot.xlabel(xColumn, size = 14)
		matplotlib.pyplot.ylabel(yColumn, size = 14)
		if 'JD' in xColumn:
			JDoffset = int(x_values[0])
			x_values = [x - JDoffset for x in x_values]
			matplotlib.pyplot.xlabel(xColumn + " - " + str(JDoffset), size = 14)
			matplotlib.pyplot.xlim(0.7, 0.9)
			matplotlib.pyplot.ylim(0.0, 4.0)

		matplotlib.pyplot.errorbar(x_values, y_values, color=c, yerr=y_errors, fmt = '.', ecolor='lightgrey', capsize=0)

	if arg.invert:
		matplotlib.pyplot.gca().invert_yaxis()
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block = False)
	fig.savefig('lightcurves.eps',dpi=100, format='eps')
	fig.savefig('lightcurves.png')
