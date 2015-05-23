#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import scipy.optimize

def function(x):
	y = a1 / (1 + math.exp(-a2*(x-a3))) + a4 + a5 * (x - a3)
	return y
	
def calcChiSquared(xValues, yValues):
	cs = 0
	for x, y in zip(xValues,yValues):
		yFit = function(x)
		print x, y, yFit
		cs+= (y-yFit)**2
	return cs
	
def func(x, a1, a2, a3, a4, a5):
	y = a1 / (1 + numpy.exp(-a2*(x-a3))) + a4 + a5 * (x - a3)
	return y
	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Fits a sigmoid function to an eclipse ingress.')
	parser.add_argument('inputfile', type=str, help='Input data in CSV format')
	 
	arg = parser.parse_args()
	print arg
	c = 'k'
	filename = arg.inputfile
	columnNames, photometry = loadingSavingUtils.loadNewCSV(filename)
	""" Data is now loaded 
	"""
	
	# Fit parameters
	a1 = -2
	a2 = 30000
	a3 = 0.81037
	a4 = 2
	a5 = -0.03

	
	matplotlib.pyplot.figure(figsize=(12, 5))
	xInput = numpy.arange(0.805, 0.815, .00001)
	y = [function(x) for x in xInput]
	matplotlib.pyplot.plot(xInput, y, color = c)
	matplotlib.pyplot.show(block= False)
	
	
	xColumn = columnNames[0]
	yColumn = columnNames[1]
	yErrors = columnNames[2]
	
	matplotlib.pyplot.figure(figsize=(12, 5))
	
	x_values = photometry[xColumn]
	y_values = photometry[yColumn]
	y_errors = photometry[yErrors]
	
	JDoffset = int(x_values[0])
	x_values = [x - JDoffset for x in x_values]
	
	
	# Make a guess of a3
	# Make a guess of a2
	stepHeight = max(y_values) - min(y_values)
	print "Step height [a1] is", stepHeight
	a1 = -stepHeight
	print "Sharpness [a2] is", a2
	
	midX = a3
	print "Midpoint of dropin x [a3] is", midX
	a3 = midX
	topHeight = max(y_values)
	print "Top height [a4] is", topHeight
	a4 = topHeight
	print "Linear slope [a5] is", a5
	
	chiSquared = calcChiSquared(x_values, y_values)
	print "chi squared:", chiSquared
	
	aguess = numpy.array([a1, a2, a3, a4, a5])
	print "a0", aguess
	aresult = scipy.optimize.curve_fit(func, x_values, y_values, aguess, y_errors)
	parameters = aresult[0]
	
	a1 = parameters[0]
	a2 = parameters[1]
	a3 = parameters[2]
	a4 = parameters[3]
	a5 = parameters[4]
	
	y_fit = [function(x) for x in x_values]
	
	print "Eclipse ingress time:", a3, a3+JDoffset
	matplotlib.pyplot.xlabel(xColumn, size = 14)
	matplotlib.pyplot.ylabel('Relative counts', size = 14)
	
	matplotlib.pyplot.xlabel(xColumn + " - " + str(JDoffset), size = 14)
			
	matplotlib.pyplot.errorbar(x_values, y_values, color=c, yerr=y_errors, fmt = '.', ecolor=c, capsize=0)
	matplotlib.pyplot.plot(x_values, y_fit, color = 'r')
	matplotlib.pyplot.draw()
	
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('ingress.eps',dpi=100, format='eps')
	fig.savefig('ingress.png',dpi=100, format='png')
	
	
	sys.exit()
	
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
