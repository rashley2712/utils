#!/usr/bin/env python
import numpy, math
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import timeClasses
import ppgplot

def filterOutNaNs(data):
	newData = {}
	for key in data.keys():
		print key
		newData[key] = []
	for index, d in enumerate(data[key]):
		for key in data.keys():
			value = data[key][index]
			if not math.isnan(value):
				newData[key].append(value)
	return newData
		
			

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses pgplot to plot light curves from rashley''s format CSV file.')
	parser.add_argument('inputfiles', type=str, nargs='+', help='Input data in CSV format')
	parser.add_argument('--bin', type=int, default = 1, help='Binning factor')
	parser.add_argument('--zero', action = 'store_true', help='Remove the mean value from the plots.... Centering around zero.')
	parser.add_argument('--errors', action = 'store_true', help='Load and plot the error bars.')
	parser.add_argument('--phaseplot', action = 'store_true', help = 'Do a phased plot')
	parser.add_argument('-e', '--ephemeris', type=str, help='Optional ephemeris file')
	 
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
	
	for filename in arg.inputfiles:
		columnNames, photometry = loadingSavingUtils.loadNewCSV(filename)
		allData.append(photometry)
	
	""" Data is now loaded 
	"""
	
	sizePerPlot = 4

	xColumn = columnNames[0]
	yColumn = columnNames[1]
	yErrors = columnNames[2]
	
	
	xLabel = xColumn
	yLabel = "flux ratio"
	
	# Initialise the plot environment 
	plotDevices = ["/xs", "pgplot.eps/ps"]
	for plotDevice in plotDevices:
		mainPGPlotWindow = ppgplot.pgopen(plotDevice)	
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pgask(True)
		
		for index, photometry in enumerate(allData):
			photometry = filterOutNaNs(photometry)
			x_values = photometry[xColumn]
			y_values = photometry[yColumn]
			y_errors = photometry[yErrors]
			x_lower, x_upper = (min(x_values), max(x_values))
			numpoints = len(x_values)
			x_range = x_upper - x_lower
			x_spacing = x_range / numpoints
			new_spacing = x_range * 86400.
			print "x-range:", x_range, " spacing:", x_spacing
			x_offset_values = [new_spacing * (x-x_lower) for x in x_values]
			ppgplot.pgsci(1)
			ppgplot.pgenv(min(x_offset_values), max(x_offset_values), min(y_values), max(y_values), 0, 0)
			ppgplot.pgslw(7)
			ppgplot.pgpt(x_offset_values, y_values, 1)
			ppgplot.pgslw(1)
			ppgplot.pgerrb(2, x_offset_values, y_values, y_errors, 0)
			ppgplot.pgerrb(4, x_offset_values, y_values, y_errors, 0)
			ppgplot.pglab(xColumn, yColumn, "")
			
