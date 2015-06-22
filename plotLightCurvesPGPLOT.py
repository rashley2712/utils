#!/usr/bin/env python
import numpy, math
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import timeClasses
import ppgplot
		

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
	
	
	
	
	c = 'k'
	sizePerPlot = 4

	xColumn = columnNames[0]
	yColumn = columnNames[1]
	yErrors = columnNames[2]
	
	xLabel = xColumn
	yLabel = "flux ratio"
	
	# Initialise the plot environment 
	plotDevices = ["/xs", "pgplot.eps/ps"]
	mainPGPlotWindow = ppgplot.pgopen(plotDevices[0])	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgask(True)
	ppgplot.pglab(xColumn, yColumn, "")
	
	for index, photometry in enumerate(allData):
		subPlot = index+1
		x_values = photometry[xColumn]
		y_values = photometry[yColumn]
		y_errors = photometry[yErrors]
		axes = matplotlib.pyplot.subplot(numPlots, 1, subPlot)
		print "subplot", numPlots, 1, subPlot
		
		matplotlib.pyplot.xlabel(xColumn, size = 14)
		matplotlib.pyplot.ylabel('Flux ratio', size = 14)
		if 'JD' in xColumn:
			JDoffset = int(x_values[0])
			x_values = [x - JDoffset for x in x_values]
			matplotlib.pyplot.xlabel(xColumn + " - " + str(JDoffset), size = 14)
			matplotlib.pyplot.xlim(0.7, 0.9)
			matplotlib.pyplot.ylim(0.0, 4.0)
			
		matplotlib.pyplot.errorbar(x_values, y_values, color=c, yerr=y_errors, fmt = '.', ecolor=c, capsize=0)
	
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('lightcurves.eps',dpi=100, format='eps')
	fig.savefig('lightcurves.png',dpi=100, format='png')
	
	
	"""	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
		else:
			matplotlib.pyplot.xlabel('Phase' , size = 14)
			
		matplotlib.pyplot.errorbar(x_values, y_values, color=c, yerr=y_errors, fmt = '.', ecolor=c, capsize=0)
		
		axes = matplotlib.pyplot.subplot(4, 1, 3)
		matplotlib.pyplot.errorbar(x_values, comparisony_values, color=c, yerr=comparisony_errors, fmt = '.', ecolor=c, capsize=0)
		
		axes = matplotlib.pyplot.subplot(4, 1, 2)
		matplotlib.pyplot.errorbar(x_values, fluxRatioValues, color=c, yerr=fluxRatioErrors, fmt = '.', ecolor=c, capsize=0)
		
		axes = matplotlib.pyplot.subplot(4, 1, 1)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.errorbar(x_values, magnitudeValues, color=c, yerr=magnitudeErrors, fmt = '.', ecolor=c, capsize=0)
		
		# Store the data for later plots
		plotData[c + 'Time'] = x_values
		plotData[c + 'Magnitudes'] = magnitudeValues
		plotData[c + 'MagnitudeErrors'] = magnitudeErrors
		plotData[c + 'FluxRatios'] = fluxRatioValues
		plotData[c + 'FluxRatioErrors'] = fluxRatioErrors
		
		matplotlib.pyplot.show()
		
	
	height = len(colours) * 4 + 1
	matplotlib.pyplot.figure(figsize=(12, height))
	
	
	for index, c in enumerate(colours):
		axes = matplotlib.pyplot.subplot(len(colours), 1, len(colours) - index)
		
		x_values = plotData[c + 'Time']
		if arg.m:
			y_values = plotData[c + 'Magnitudes']
			y_errors = plotData[c + 'MagnitudeErrors']
		else: 
			y_values = plotData[c + 'FluxRatios']
			y_errors = plotData[c + 'FluxRatioErrors']
			
		
		if arg.m: 
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.ylabel(r"$%s_{mag}$"%filters[c], size = 18)
		else:
			matplotlib.pyplot.ylabel(r"$%s_{flux}$"%filters[c], size = 18)
			 
		matplotlib.pyplot.errorbar(x_values, y_values, color=c, yerr=y_errors, fmt = '.', ecolor=c, capsize=0)
		

		if index==0: 
			if (not arg.phaseplot):
				matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
			else:
				matplotlib.pyplot.xlabel('Phase' , size = 14)

	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('lightcurves.eps',dpi=100, format='eps')
	fig.savefig('lightcurves.png',dpi=100, format='png')

	for c in colours:

		filename = arg.outcsv + "_" + c + ".csv"
		data = []
		times = plotData[c + 'Time']
		fluxRatioValues = plotData[c + 'FluxRatios']
		fluxRatioErrors = plotData[c + 'FluxRatioErrors']
		for index, time in enumerate(times):
			record = {}
			record['MJD'] = time
			record['fluxRatio'] = fluxRatioValues[index]
			record['fluxRatioError'] = fluxRatioErrors[index]
			data.append(record)
		loadingSavingUtils.writeSingleChannelCSV(filename, data)
		
		
	if len(arg.channels)==1:
		sys.exit()
	
	numPlots = len(arg.colourplots)
	
	matplotlib.pyplot.figure(figsize=(12, numPlots*4 + 1))
	print "numplots", numPlots
		
	for colourPlotRequested in arg.colourplots:
		print "Colour plot requested", colourPlotRequested
		if colourPlotRequested=='gr':
			axes = matplotlib.pyplot.subplot(numPlots, 1, numPlots)
			rx_values = plotData['rTime']
			r_values = plotData['rMagnitudes']
			r_errors = plotData['rMagnitudeErrors']
			gx_values = plotData['gTime']
			g_values = plotData['gMagnitudes']
			g_errors = plotData['gMagnitudeErrors']
			
			gMinusrValues = []
			gMinusrErrors = []
			x_values = []
			for index, r in enumerate(rx_values):
				time = r
				try:
					gIndex = gx_values.index(time)
				except ValueError:
					print "Couldn't find a corresponding data point in 'g' for the 'r' data at", time
					continue
					
				g_value = g_values[gIndex]
				gMinusr = g_value - r_values[index]
				gMinusrError = math.sqrt( r_errors[index]**2 + g_errors[gIndex]**2 )
				x_values.append(time)
				gMinusrValues.append(gMinusr)
				gMinusrErrors.append(gMinusrError)
				
			
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.errorbar(x_values, gMinusrValues, color='k', yerr=gMinusrErrors, fmt = '.', ecolor='k', capsize=0)
			matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
		
			matplotlib.pyplot.ylabel(r"$(g-i)_{mag}$", size = 18)
			
			
		if colourPlotRequested=='ug':
			axes = matplotlib.pyplot.subplot(numPlots, 1, 1)
			bx_values = plotData['bTime']
			b_values = plotData['bMagnitudes']
			b_errors = plotData['bMagnitudeErrors']
			gx_values = plotData['gTime']
			g_values = plotData['gMagnitudes']
			g_errors = plotData['gMagnitudeErrors']
			
			uMinusgValues = []
			uMinusgErrors = []
			x_values = []
			for index, b in enumerate(bx_values):
				
				time = b
				gIndex, distance = statsUtils.findNearestTime(time, gx_values, g_values)
				print "Closest g-time", gx_values[gIndex], "to b-time", time, distance
					
				g_value = g_values[gIndex]
				b_value = b_values[index]
				uMinusg = b_value - g_value
				uMinusgError = math.sqrt( b_errors[index]**2 + g_errors[gIndex]**2 )
				x_values.append(time)
				uMinusgValues.append(uMinusg)
				uMinusgErrors.append(uMinusgError)
				#print bx_values[bIndex], b_value, time, g_values[index] 
				
			
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.errorbar(x_values, uMinusgValues, color='k', yerr=uMinusgErrors, fmt = '.', ecolor='k', capsize=0)
			matplotlib.pyplot.ylabel(r"$(u-g)_{mag}$", size = 18)
	
			
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('colourcurves.eps',dpi=100, format='eps')
	fig.savefig('colourcurves.png',dpi=100, format='png')
	"""
