#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import scipy.optimize
import time
import random
import ppgplot

def function(x):
	y = a1 / (1 + numpy.exp(-a2*(x-a3))) + a4 + a5 * (x - a3)
	return y
	
def calcChiSquared(xValues, yValues):
	cs = 0
	for x, y in zip(xValues,yValues):
		yFit = function(x)
		cs+= (y-yFit)**2
	return cs
	
def func(x, a1, a2, a3, a4, a5):
	y = a1 / (1 + numpy.exp(-a2*(x-a3))) + a4 + a5 * (x - a3)
	return y
	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Fits a sigmoid function to an eclipse ingress.')
	parser.add_argument('inputfile', type=str, help='Input data in CSV format')
	parser.add_argument('--trim', type = int, default=5, help='Max size of trimming of the data set. [default: 5].') 
	parser.add_argument('-n', '--iterations', type = int, default=1000, help='Number of iterations for the MC bootstrap method. [default: 1000].') 
	arg = parser.parse_args()
	print arg
	runStr = arg.inputfile[:6]
	trimSize = arg.trim
	c = 'k'
	filename = arg.inputfile
	columnNames, photometry = loadingSavingUtils.loadNewCSV(filename)
	""" Data is now loaded 
	"""
	xColumn = columnNames[0]
	yColumn = columnNames[1]
	yErrors = columnNames[2]
	x_values = photometry[xColumn]
	y_values = photometry[yColumn]
	y_errors = photometry[yErrors]
	JDoffset = int(x_values[0])
	x_values = [x - JDoffset for x in x_values]
	
	# Determine ingress or egress
	beginY = numpy.mean(y_values[0:4])
	endY = numpy.mean(y_values[-5:-1])
	print "begin y", beginY, "end y", endY
	if endY>beginY:
		egress = True
		print "This is an egress."
		fileappendix = "egress"
	else:
		egress = False
		print "This is an ingress."
		fileappendix = "ingress"
		
	# Initial parameters
	drop = endY - beginY
	a1 = drop
	print "Drop is [a0]:", a1
	a2 = 30000.
	print "Sharpness [a2] is", a2
	a3 = numpy.median(x_values)
	print "Midpoint of drop is [a3]", a3
	a4 = beginY
	print "Initial value is [a4]", a4
	a5 = 0.
	print "Linear slope [a5] is", a5
	
	# Do a test plot first
	y_fit = function(x_values)
	matplotlib.pyplot.figure(figsize=(8, 5))
	matplotlib.pyplot.plot(x_values, y_fit)
	matplotlib.pyplot.plot(x_values, y_values)
	matplotlib.pyplot.show(block=False)
	matplotlib.pyplot.figure(figsize=(8, 5))
	

	
	chiSquared = calcChiSquared(x_values, y_values)
	#print "chi squared:", chiSquared
	
	aguess = numpy.array([a1, a2, a3, a4, a5])
	startGuess = aguess
	print "a0", aguess
	aresult = scipy.optimize.curve_fit(func, x_values, y_values, aguess, y_errors)
	parameters = aresult[0]
	originalResult = aresult[0]
	
	a1 = parameters[0]
	a2 = parameters[1]
	a3 = parameters[2]
	a4 = parameters[3]
	a5 = parameters[4]
	
	y_fit = [function(x) for x in x_values]
	lowerX = min(x_values)
	upperX = max(x_values)
	
	print "Eclipse %s time: %f or %5.8f"%(fileappendix, a3, a3+JDoffset)
	
	matplotlib.pyplot.xlabel(xColumn, size = 14)
	matplotlib.pyplot.ylabel('Relative counts', size = 14)
	
	matplotlib.pyplot.xlabel(xColumn + " - " + str(JDoffset), size = 14)
			
	matplotlib.pyplot.errorbar(x_values, y_values, color=c, yerr=y_errors, fmt = '.', ecolor=c, capsize=0)
	#matplotlib.pyplot.plot(x_values, y_fit, color = 'r')
	yLims = matplotlib.pyplot.gca().get_ylim()
	xLims = matplotlib.pyplot.gca().get_xlim()
	print "x, y-limits: ", xLims, yLims
	steps = 1000
	size = (xLims[1]-xLims[0])/steps
	xFit = numpy.arange(xLims[0], xLims[1], size)
	yFit = [function(x) for x in xFit]
	matplotlib.pyplot.plot(xFit, yFit, color = 'g')
	
	matplotlib.pyplot.plot([a3, a3], [yLims[0], yLims[1]], color = 'g', linestyle='dashed')
	
	fig = matplotlib.pyplot.gcf()
	fig.suptitle(runStr + ' ' + fileappendix, fontsize=20)

	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block = False)
	fig.savefig(runStr +'_%s.eps'%fileappendix,dpi=100, format='eps')
	fig.savefig(runStr +'_%s.png'%fileappendix,dpi=100, format='png')
	
	plotDevices = ["/xs", "%s.eps/ps"%fileappendix]
	for plotDevice in plotDevices:
		mainPGPlotWindow = ppgplot.pgopen(plotDevice)
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pgpap(10, 0.618)
		ppgplot.pgsci(1)
		ppgplot.pgenv(min(x_values), max(x_values), yLims[0], yLims[1], 0, 0)
		ppgplot.pgslw(7)
		ppgplot.pgpt(x_values, y_values, 1)
		ppgplot.pgslw(1)
		ppgplot.pgerrb(2, x_values, y_values, y_errors, 0)
		ppgplot.pgerrb(4, x_values, y_values, y_errors, 0)
		ppgplot.pgsls(2)
		ppgplot.pgline(xFit, yFit)
		ppgplot.pgsls(3)
		ppgplot.pgline([a3, a3], [yLims[0], yLims[1]])
		ppgplot.pgsls(1)
		ppgplot.pglab(xColumn + " - " + str(JDoffset), "flux ratio", "")
		ppgplot.pgclos()
	time.sleep(3)
	
	times = []
	sharpness = []
	random.seed()
	for n in range(arg.iterations):
		# Give all the points a random bump...
		y_perturbed = []
		for y, y_error in zip(y_values, y_errors):
			y_p = numpy.random.normal(y, y_error)
			y_perturbed.append(y_p)	
		y_perturbed = numpy.array(y_perturbed)
		# Trim off a few points from each end of the data. 
		leftTrim = random.randrange(0, trimSize)
		rightTrim = random.randrange(0, trimSize)
		trimmed_x = x_values[leftTrim: len(x_values)-rightTrim]
		trimmed_y = y_values[leftTrim: len(x_values)-rightTrim]
		trimmed_yp = y_perturbed[leftTrim: len(x_values)-rightTrim]
		trimmed_errors = y_errors[leftTrim: len(x_values)-rightTrim]
		# Run a new fit...
		aguess = numpy.array([a1, a2, a3, a4, a5])
		aguess = startGuess
		aresult = scipy.optimize.curve_fit(func, trimmed_x, trimmed_yp, aguess, trimmed_errors)
		parameters = aresult[0]
		a1 = parameters[0]
		a2 = parameters[1]
		a3 = parameters[2]
		a4 = parameters[3]
		a5 = parameters[4]
		print "Montecarlo test number:", n, "Trim: [%d, %d]"%(leftTrim, rightTrim),"Eclipse ingress time: %f or %5.10f"%(a3, a3+JDoffset)
		times.append(a3)
		sharpness.append(a2)
		
		steps = 1000
		size = (xLims[1]-xLims[0])/steps
		x_inter = numpy.arange(xLims[0], xLims[1], size)
		y_fit = [function(x) for x in x_inter]
		
		matplotlib.pyplot.clf()
		matplotlib.pyplot.xlabel(xColumn, size = 14)
		matplotlib.pyplot.ylabel('Relative counts', size = 14)
		matplotlib.pyplot.xlabel(xColumn + " - " + str(JDoffset), size = 14)
		matplotlib.pyplot.errorbar(trimmed_x, trimmed_y, color=c, yerr=trimmed_errors, fmt = '.', ecolor=c, capsize=0)
		matplotlib.pyplot.scatter(trimmed_x, trimmed_yp, color='r', marker = 'o')
		matplotlib.pyplot.plot([a3, a3], [yLims[0], yLims[1]], color = 'g', linestyle='dashed')
		matplotlib.pyplot.plot(x_inter, y_fit, color = 'r')
		matplotlib.pyplot.ylim(yLims)
		matplotlib.pyplot.xlim(xLims)
		matplotlib.pyplot.draw()
		#time.sleep(0.2)
		
	print "Time [a3]:"
	print "Mean: %5.10f  Stddev:%2.12f or %f seconds"%(numpy.mean(times), numpy.std(times), numpy.std(times)*86400.)	
	print "Original time: %5.8f, Montecarlo mean: %5.8f, stddev: %2.12f"%(originalResult[2] + JDoffset, numpy.mean(times) + JDoffset, numpy.std(times))
	print "Sharpness [a2]:"
	print "Mean: %5.10f  Stddev:%2.12f"%(numpy.mean(sharpness), numpy.std(sharpness))	
	fig = matplotlib.pyplot.gcf()
	#matplotlib.pyplot.show()
	"""
	fig.savefig('ingress.eps',dpi=100, format='eps')
	fig.savefig('ingress.png',dpi=100, format='png')
	"""
	
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
