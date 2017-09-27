#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import timeClasses
from matplotlib import rc
		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses matplotlib to plot light curves from Liverpool telescope data.')
	parser.add_argument('inputfiles', type=str, nargs='+', help='Input data in TSV format')
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
	
	numPlots = 0
	for filename in arg.inputfiles:
		columnNames, photometry = loadingSavingUtils.loadLiverpoolTSV(filename)
		allData.append(photometry)
		numPlots+= 1
	
	
	
	""" Data is now loaded 
	"""
	c = 'k'
	sizePerPlot = 6

	for photometry in allData:
		# Add the heliocentric correction to the times
		JDs = photometry['JD']
		corrections = photometry['HC']
		HJDs = [ jd + corr for jd, corr in zip(JDs, corrections)]
		photometry['HJD'] = HJDs
		
	if arg.bin!=1:
		binsize = arg.bin
		for photometry in allData:
			HJDs = photometry['HJD']
			flux = photometry['dmag']
			err = photometry['err']
			for d, f, e in zip(HJDs, flux, err):
				print d, f, e
			startIndex = 0
			endIndex = len(HJDs) - binsize
			print startIndex, endIndex
		
		binnedHJDs = []
		binnedFluxes = []
		binnedErrors = []
			
		for index in range(startIndex, endIndex+binsize, binsize):
			dates = 0
			fluxes = 0
			errors = 0
			for step in range(binsize):
				dates+=HJDs[index+step]
				fluxes+=flux[index+step]
				errors+=err[index+step]
			dates/=binsize
			fluxes/=binsize
			errors/=binsize
			binnedHJDs.append(dates)
			binnedFluxes.append(fluxes)
			binnedErrors.append(errors)
		photometry['HJD'] = binnedHJDs
		photometry['dmag'] = binnedFluxes
		photometry['err'] = binnedErrors	
	
	#sys.exit()	
	
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	## for Palatino and other serif fonts use:
	rc('font',**{'family':'serif','serif':['Palatino']})
	rc('text', usetex=True)

	xColumn = 'HJD'
	yColumn = 'dmag'
	yErrors = 'err'
	figSize = 10
	matplotlib.pyplot.figure(figsize=(figSize, figSize / 1.618))
	
	labelSize = 20
	tickSize = 18
	for index, photometry in enumerate(allData):
		subPlot = index+1
		x_values = photometry[xColumn]
		y_values = photometry[yColumn]
		y_errors = photometry[yErrors]
		axes = matplotlib.pyplot.subplot(numPlots, 1, subPlot)
		print "subplot", numPlots, 1, subPlot
		
		matplotlib.pyplot.xlabel(xColumn, size = labelSize)
		matplotlib.pyplot.ylabel('Apparent magnitude', size = labelSize)
		if 'JD' in xColumn:
			JDoffset = int(x_values[0])
			x_values = [x - JDoffset for x in x_values]
			matplotlib.pyplot.xlabel(xColumn + " - " + str(JDoffset), size = labelSize)
			#matplotlib.pyplot.xlim(0.7, 0.9)
			#matplotlib.pyplot.ylim(0.0, 4.0)
			
		matplotlib.pyplot.errorbar(x_values, y_values, color=c, yerr=y_errors, fmt = '.', ecolor='0.75', capsize=0)
		
		for label in (axes.get_xticklabels() + axes.get_yticklabels()):
			#label.set_fontname('Arial')
			label.set_fontsize(tickSize)
		matplotlib.pyplot.gca().invert_yaxis()
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	# fig.savefig('lightcurves.eps',dpi=100, format='eps')
	fig.savefig('lightcurves.pdf',dpi=100, format='pdf')
	
	
