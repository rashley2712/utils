#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils

		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Converts MJDs to Phase given a certain ephemeris.')
	parser.add_argument('inputdata', type=str, help='Input data in CSV format')
	parser.add_argument('--channels', nargs='+', default = ['r', 'g', 'b'], type=str, help='Channels to process')
	parser.add_argument('--outcsv', type=str, default = "temp_phased.csv", help='Output each channel to csv with this name.')
	parser.add_argument('--T0', type=float, default = 0, help='T0 time of the ephemeris.')
	parser.add_argument('--E1', type=float, default = 0, help='E1 period of the ephemeris.')
	arg = parser.parse_args()
	print arg

	coloursRequested = arg.channels

	if (arg.T0==0):
		t0 = 55126.3960
	else:
		t0 = arg.T0

	if (arg.E1==0):
		e1 = 0.081376
	else:
		e1 = arg.E1
		
	 
	arg = parser.parse_args()
	print arg
	
	photometry = loadingSavingUtils.loadCSV(arg.inputdata)
	
	filters = { 'r':'i', 'g':'g', 'b':'u' }
	
	""" Data is now loaded 
	"""

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
		data = photometry[c]
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
			#print fluxRatio, fluxRatioError
			magnitude = -2.5 * math.log10(fluxRatio)
			magnitudeError = -2.5 * fluxRatioError / math.log(10) / fluxRatio
			magnitudeValues.append(magnitude)
			magnitudeErrors.append(magnitudeError)
			
			duration = d['MJD'] - t0
			epochs = math.floor(duration/e1)
			phase = (duration - (epochs * e1))/e1
			phase2 = duration/e1 - epochs 
			print d['MJD'], duration, epochs, phase, phase2
			d['Phase'] = phase2
			d['Epoch'] = epochs
			#print fluxRatio, fluxRatioError, magnitude, magnitudeError
			
		
		
		magnitudeMean = numpy.mean(magnitudeValues)
		magnitudeValues = [m - magnitudeMean for m in magnitudeValues]
		
		outputFilename = arg.outcsv
		loadingSavingUtils.writeCSV(outputFilename, photometry, colours = coloursRequested)	
		
		