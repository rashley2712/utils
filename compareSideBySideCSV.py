#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils

		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses matplotlib to plot light curves from rashley''s format CSV file.')
	parser.add_argument('inputdata1', type=str, help='Input data in CSV format')
	parser.add_argument('inputdata2', type=str, help='Input data in CSV format')
	parser.add_argument('-c', nargs='+', default = ['Counts_1'], type=str, help='Columns to plot')
	parser.add_argument('--channels', nargs='+', default = ['r', 'g', 'b'], type=str, help='Channels to plot')
	parser.add_argument('-m', action = 'store_true', help='Use magnitude scale')
	parser.add_argument('--bin', type=int, default = 1, help='Binning factor')
	parser.add_argument('--outcsv', type=str, default = "temp", help='Output each channel to csv with this name.')
	parser.add_argument('--zero', action = 'store_true', help='Remove the mean value from the plots.... Centering around zero.')
	parser.add_argument('--errors', action = 'store_true', help='Load and plot the error bars.')
	parser.add_argument('--colourplots', nargs='+', default = ['gr'], type=str, help='Colour plots requested')
	
	 
	arg = parser.parse_args()
	print arg
	
	photometry1 = loadingSavingUtils.loadCSV(arg.inputdata1)
	photometry2 = loadingSavingUtils.loadCSV(arg.inputdata2)
	
	
	""" Data is now loaded 
	"""
	filters = { 'r':'i', 'g':'g', 'b':'u' }
	
	plotData = {}
	colours = arg.channels
	for c in colours:
		x_values = []
		y_values = []
		y_errors = []
		comparisony_values = []
		comparisony_errors = []
		fluxRatioValues = []
		fluxRatioErrors = []
		magnitudeValues = []
		magnitudeErrors = []
		data = photometry1[c]
		for d in data:
			x_values.append(d['MJD'])
			y_values.append(d['Counts_1'])
			y_errors.append(d['Sigma_1'])
			comparisony_values.append(d['Counts_2'])
			comparisony_errors.append(d['Sigma_2'])
			fluxRatio = d['Counts_1'] / d['Counts_2']
			fluxRatioError = fluxRatio * math.sqrt( (d['Sigma_1']/d['Counts_1'])**2 + (d['Sigma_2']/d['Counts_2'])**2 )
			fluxRatioValues.append(fluxRatio)
			fluxRatioErrors.append(fluxRatioError)
			magnitude = -2.5 * math.log10(fluxRatio)
			magnitudeError = -2.5 * fluxRatioError / math.log(10) / fluxRatio
			magnitudeValues.append(magnitude)
			magnitudeErrors.append(magnitudeError)
			
		
		if (arg.bin!=1):
			x_temp, y_values, y_errors = statsUtils.binDataWithErrors(x_values, y_values, y_errors, arg.bin)
			x_temp, fluxRatioValues, fluxRatioErrors = statsUtils.binDataWithErrors(x_values, fluxRatioValues, fluxRatioErrors, arg.bin)
			x_temp, magnitudeValues, magnitudeErrors = statsUtils.binDataWithErrors(x_values, magnitudeValues, magnitudeErrors, arg.bin)
			x_values, comparisony_values, comparisony_errors = statsUtils.binDataWithErrors(x_values, comparisony_values, comparisony_errors, arg.bin)
		
		magnitudeMean = numpy.mean(magnitudeValues)
		magnitudeValues = [m - magnitudeMean for m in magnitudeValues]
		
		MJDoffset = int(min(x_values))
		print "MJD offset:",MJDoffset
		
		x_values = [x - MJDoffset for x in x_values]
		
		matplotlib.pyplot.figure(figsize=(12, 6))
		axes = matplotlib.pyplot.subplot(1, 1, 1)
		matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)
		
		matplotlib.pyplot.plot(x_values, fluxRatioValues,  c+'s', label = 'traditional pipeline')
		matplotlib.pyplot.ylabel(r"$Counts_{NN Ser}/Counts_{comparison}$", size = 18)
	
			
		#matplotlib.pyplot.errorbar(x_values, fluxRatioValues, color=c, yerr=fluxRatioErrors, fmt = '.', ecolor=c)
		
		# Store the data for later plots
		plotData[c + 'Time'] = x_values
		plotData[c + 'Magnitudes'] = magnitudeValues
		plotData[c + 'MagnitudeErrors'] = magnitudeErrors
		
		x_values = []
		y_values = []
		y_errors = []
		comparisony_values = []
		comparisony_errors = []
		fluxRatioValues = []
		fluxRatioErrors = []
		magnitudeValues = []
		magnitudeErrors = []
		data = photometry1[c]
		data = photometry2[c]
		for d in data:
			x_values.append(d['MJD'])
			y_values.append(d['Counts_1'])
			y_errors.append(d['Sigma_1'])
			comparisony_values.append(d['Counts_2'])
			comparisony_errors.append(d['Sigma_2'])
			fluxRatio = d['Counts_1'] / d['Counts_2']
			print d['MJD'], fluxRatio
			fluxRatioError = fluxRatio * math.sqrt( (d['Sigma_1']/d['Counts_1'])**2 + (d['Sigma_2']/d['Counts_2'])**2 )
			fluxRatioValues.append(fluxRatio)
			fluxRatioErrors.append(fluxRatioError)
			magnitude = -2.5 * math.log10(fluxRatio)
			magnitudeError = -2.5 * fluxRatioError / math.log(10) / fluxRatio
			magnitudeValues.append(magnitude)
			magnitudeErrors.append(magnitudeError)
			
		
		if (arg.bin!=1):
			x_temp, y_values, y_errors = statsUtils.binDataWithErrors(x_values, y_values, y_errors, arg.bin)
			x_temp, fluxRatioValues, fluxRatioErrors = statsUtils.binDataWithErrors(x_values, fluxRatioValues, fluxRatioErrors, arg.bin)
			x_temp, magnitudeValues, magnitudeErrors = statsUtils.binDataWithErrors(x_values, magnitudeValues, magnitudeErrors, arg.bin)
			x_values, comparisony_values, comparisony_errors = statsUtils.binDataWithErrors(x_values, comparisony_values, comparisony_errors, arg.bin)
		
		magnitudeMean = numpy.mean(magnitudeValues)
		magnitudeValues = [m - magnitudeMean for m in magnitudeValues]
		
		x_values = [x - MJDoffset for x in x_values]
		top = numpy.max(fluxRatioValues)
		bottom = numpy.min(fluxRatioValues)
		offset = (top-bottom) * 0.2
		print "Offset applied:", offset
		fluxRatioValues = [f + offset for f in fluxRatioValues]
		
		matplotlib.pyplot.plot(x_values, fluxRatioValues,  c+'o', label = 'automated pipeline')
		matplotlib.pyplot.legend(loc='lower left')
		
		matplotlib.pyplot.show()
		
	""" Now do the stacked plot
	"""
	
	height = len(colours) * 4 + 1
	matplotlib.pyplot.figure(figsize=(12, height))
	
	
	for index, c in enumerate(colours):
		axes = matplotlib.pyplot.subplot(len(colours), 1, len(colours) - index)
		
		x_values = plotData[c + 'Time']
		y_values = plotData[c + 'Magnitudes']
		y_errors = plotData[c + 'MagnitudeErrors']
		
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.errorbar(x_values, y_values, color=c, yerr=y_errors, fmt = '.', ecolor=c)
		matplotlib.pyplot.ylabel(r"$%s_{mag}$"%filters[c], size = 18)
	

		if index==0: matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset), size = 14)

	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('lightcurves.eps',dpi=100, format='eps')
	fig.savefig('lightcurves.png',dpi=100, format='png')
	
	filename = arg.outcsv + "_r.csv"
	data = []
	times = plotData['rTime']
	magnitudes = plotData['rMagnitudes']
	magnitudeErrors = plotData['rMagnitudeErrors']
	for index, time in enumerate(times):
		record = {}
		record['MJD'] = time
		record['magnitude'] = magnitudes[index]
		record['magnitudeError'] = magnitudeErrors[index]
		data.append(record)
	loadingSavingUtils.writeSingleChannelCSV(filename, data)
	
	filename = arg.outcsv + "_g.csv"
	data = []
	times = plotData['gTime']
	magnitudes = plotData['gMagnitudes']
	magnitudeErrors = plotData['gMagnitudeErrors']
	for index, time in enumerate(times):
		record = {}
		record['MJD'] = time
		record['magnitude'] = magnitudes[index]
		record['magnitudeError'] = magnitudeErrors[index]
		data.append(record)
	loadingSavingUtils.writeSingleChannelCSV(filename, data)
		
		
	""" Now do the colour plot(s)
	"""
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
			matplotlib.pyplot.errorbar(x_values, gMinusrValues, color='k', yerr=gMinusrErrors, fmt = '.', ecolor='k')
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
			matplotlib.pyplot.errorbar(x_values, uMinusgValues, color='k', yerr=uMinusgErrors, fmt = '.', ecolor='k')
			matplotlib.pyplot.ylabel(r"$(u-g)_{mag}$", size = 18)
	
			
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('colourcurves.eps',dpi=100, format='eps')
	fig.savefig('colourcurves.png',dpi=100, format='png')
	
