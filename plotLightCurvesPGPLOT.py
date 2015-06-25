#!/usr/bin/env python
import numpy, math
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import timeClasses
import ppgplot
import generalUtils

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses pgplot to plot light curves from rashley''s format CSV file.')
	parser.add_argument('inputfiles', type=str, nargs='+', help='Input data in CSV format')
	parser.add_argument('--bin', type=int, default = 1, help='Binning factor')
	parser.add_argument('--zero', action = 'store_true', help='Remove the mean value from the plots.... Centering around zero.')
	parser.add_argument('--errors', action = 'store_true', help='Load and plot the error bars.')
	parser.add_argument('--ask', action = 'store_true', help='PGPLOT will ask the user before proceeding to the next plot.')
	parser.add_argument('--phaseplot', action = 'store_true', help = 'Do a phased plot')
	parser.add_argument('-e', '--ephemeris', type=str, help='Optional ephemeris file')
	parser.add_argument('--list', action='store_true', help='Specify this option if the input file is actually a list of input files.')
	
	 
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
		
	filenames = []
	if arg.list:
		# Load the list of files.
		if len(arg.inputfiles)>1: 
			print "You can only give me one list of filenames."
			sys.exit()
		filename = arg.inputfiles[0]
		fileList = open(filename, 'r')
		for line in fileList:
			filenames.append(str(line.strip()))
	else:
		filenames = arg.inputfiles
	
	
	allData = []
	
	for filename in filenames:
		columnNames, photometry = loadingSavingUtils.loadNewCSV(filename)
		photometry = generalUtils.filterOutNaNs(photometry)
		photometry['runName'] = filename
		allData.append(photometry)

	xColumn = columnNames[0]
	yColumn = columnNames[1]
	yErrors = columnNames[2]
	
	""" Data is now loaded 
	"""
	# Sort the data by mjd
	sortedData = sorted(allData, key=lambda object: object[xColumn][0], reverse = False)
	allData = sortedData
	
	
	# Compute phases
	if hasEphemeris:
		for photometry in allData:
			phases = []
			for mjd in photometry[xColumn]:
				jd = mjd + 2400000.5
				p = ephemeris.getPhase(jd)
				if p<0.5:
					p = p + 1.0
				phases.append(p)
			photometry["phase"] = phases
	
	
	sizePerPlot = 4

	if hasEphemeris: xColumn = "phase"
	
	xLabel = xColumn
	yLabel = "flux ratio"
	
	# Find the best y-limits
	lowerY = 1E8
	upperY = -1E8
	for photometry in allData:
		yData = photometry[yColumn]
		if max(yData)>upperY: upperY = max(yData)
		if min(yData)<lowerY: lowerY = min(yData)
	
	# Initialise the plot environment 
	plotDevices = ["/xs"]
	for plotDevice in plotDevices:
		mainPGPlotWindow = ppgplot.pgopen(plotDevice)
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pgpap(10, 0.618)
		ppgplot.pgask(arg.ask)
		
		for index, photometry in enumerate(allData):
			x_values = photometry[xColumn]
			y_values = photometry[yColumn]
			y_errors = photometry[yErrors]
			x_lower, x_upper = (min(x_values), max(x_values))
			numpoints = len(x_values)
			if "JD" in xColumn:
				x_offset = int(x_lower)
				x_values = [(x-x_offset) for x in x_values]
				xLabel= xColumn + " - %d"%x_lower
			ppgplot.pgsci(1)
			ppgplot.pgenv(min(x_values), max(x_values), lowerY, upperY, 0, 0)
			ppgplot.pgslw(7)
			ppgplot.pgpt(x_values, y_values, 1)
			ppgplot.pgslw(1)
			ppgplot.pgerrb(2, x_values, y_values, y_errors, 0)
			ppgplot.pgerrb(4, x_values, y_values, y_errors, 0)
			ppgplot.pglab(xLabel, yLabel, photometry["runName"])
			
		ppgplot.pgclos()
	if not hasEphemeris:
		sys.exit()
	
				
	# Restrict the light-curve to a subset of phase
	phaseLimits = (0.95, 1.05)
	for photometry in allData:
		pnew = []
		ynew = []
		yenew = []
		for (p, y, ye) in zip(photometry["phase"], photometry[yColumn], photometry[yErrors]):
			if p>phaseLimits[0] and p<phaseLimits[1]:
				pnew.append(p)
				ynew.append(y)
				yenew.append(ye)
		photometry["phase"] = pnew
		photometry[yColumn] = ynew
		photometry[yErrors] = yenew
	
	# Find the best y-limits
	lowerY = 1E8
	upperY = -1E8
	for photometry in allData:
		yData = photometry[yColumn]
		if max(yData)>upperY: upperY = max(yData)
		if min(yData)<lowerY: lowerY = min(yData)
	
	# Initialise the plot environment 
	plotDevices = ["/xs", "eclipses.eps/ps"]
	for plotDevice in plotDevices:
		mainPGPlotWindow = ppgplot.pgopen(plotDevice)	
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pgask(arg.ask)
		
		for index, photometry in enumerate(allData):
			x_values = photometry[xColumn]
			y_values = photometry[yColumn]
			y_errors = photometry[yErrors]
			x_lower, x_upper = (min(x_values), max(x_values))
			numpoints = len(x_values)
			if "JD" in xColumn:
				x_offset = int(x_lower)
				x_values = [(x-x_offset) for x in x_values]
				xLabel+= " - %d"%x_lower
			ppgplot.pgsci(1)
			ppgplot.pgenv(phaseLimits[0], phaseLimits[1], lowerY, upperY, 0, 0)
			ppgplot.pgslw(7)
			ppgplot.pgpt(x_values, y_values, 1)
			ppgplot.pgslw(1)
			ppgplot.pgerrb(2, x_values, y_values, y_errors, 0)
			ppgplot.pgerrb(4, x_values, y_values, y_errors, 0)
			ppgplot.pglab(xLabel, yLabel, photometry["runName"])
		ppgplot.pgclos()
		
	# Plot the stacked image
	numPlots = len(allData)
	offset = 4
	upperY = offset * (numPlots - 1) + upperY
	plotDevices = ["/xs", "stacked_eclipses.eps/vps"]
	for plotDevice in plotDevices:
		mainPGPlotWindow = ppgplot.pgopen(plotDevice)	
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pgpap(6.18, 1.618)
		ppgplot.pgask(arg.ask)
		ppgplot.pgsci(1)
		ppgplot.pgenv(phaseLimits[0], phaseLimits[1], lowerY, upperY, 0, 0)
		ppgplot.pglab(xLabel, yLabel, "")
		for index, photometry in enumerate(allData):
			mjdInt = int(photometry['BMJD'][0])
			x_values = photometry[xColumn]
			y_values = photometry[yColumn]
			y_values = [y + offset*index for y in y_values]
			y_errors = photometry[yErrors]
			x_lower, x_upper = (min(x_values), max(x_values))
			numpoints = len(x_values)
			ppgplot.pgslw(7)
			ppgplot.pgpt(x_values, y_values, 1)
			ppgplot.pgslw(1)
			ppgplot.pgerrb(2, x_values, y_values, y_errors, 0)
			ppgplot.pgerrb(4, x_values, y_values, y_errors, 0)
			ppgplot.pgptxt(0.98, offset*index + offset/3, 0, 0, "MJD: %d"%mjdInt)
		ppgplot.pgclos()