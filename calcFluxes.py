#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import numpy, math
import matplotlib.pyplot
import argparse
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot

def quad(x, a1, a2, a3):
	y = a1 * x * x + a2 *x + a3
	return y
	
def linear(x, m, c):
	y = m * x + c
	return y

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a series of spectra from the JSON files. Calculates some flux measurements.')
	parser.add_argument('inputFiles', type=str, nargs='+', help='Molly files containing the spectra')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	parser.add_argument('--list', action='store_true', help='Specify this option if the input file is actually a list of input files.')
	parser.add_argument('--model', type=str, help='Optional model spectrum (to subtract).')
	 
	arg = parser.parse_args()
	print arg
	
	print "Astropy version:", astropy.__version__
	
	if arg.e!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.e)
		print ephemeris
	else:
		hasEphemeris = False

	spectra = []
	
	filenames = []
	if arg.list:
		# Load the list of files.
		if len(arg.inputFiles)>1: 
			print "You can only give me one list of filenames."
			sys.exit()
		filename = arg.inputFiles[0]
		fileList = open(filename, 'r')
		for line in fileList:
			filenames.append(str(line.strip()))
	else:
		filenames = arg.inputFiles
	print "Files to be loaded", filenames
	
	
	spectra = []
	for fileIndex, f in enumerate(filenames):
		spectrum = spectrumClasses.spectrumObject()
		spectrum.loadFromJSON(f)
		print "Loaded %s, contains %s."%(f, spectrum.objectName)
		spectra.append(spectrum)
		
	numSpectra = len(spectra)


	mainPGPlotWindow = ppgplot.pgopen('/xs')	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	yUpper = 2.5
	yLower = 0.0
	
	for spectrum in spectra:
		ppgplot.pgsci(1)
		lowerWavelength = min(spectrum.wavelengths)
		upperWavelength = max(spectrum.wavelengths)
		lowerFlux = min(spectrum.flux)
		upperFlux = max(spectrum.flux)
		ppgplot.pgenv(lowerWavelength, upperWavelength, yLower, yUpper, 0, 0)
		ppgplot.pgline(spectrum.wavelengths, spectrum.flux)
		if hasEphemeris:
			ppgplot.pglab("wavelength", "flux", "%s [%f]"%(spectrum.objectName, ephemeris.getPhase(spectrum.HJD)))
		else:
			ppgplot.pglab("wavelength", "flux", spectrum.objectName)
		
			

		# Get a sub spectrum for line fitting
		subSpectrum = copy.deepcopy(spectrum)
		subSpectrum.snipWavelengthRange(min(subSpectrum.wavelengths), 7450)
		subSpectrum.snipWavelengthRange(7550, 8130)
		subSpectrum.snipWavelengthRange(8170, 8222)
		subSpectrum.snipWavelengthRange(8262, max(subSpectrum.wavelengths))
		ppgplot.pgsci(2)
		ppgplot.pgline(subSpectrum.wavelengths, subSpectrum.flux)

		# Fit and plot a line through these three regions of the spectrum
		m = 0.0
		c = 0.0
		guess = numpy.array([m, c])
		x_values = subSpectrum.wavelengths
		y_values = subSpectrum.flux
		y_errors = numpy.ones(len(x_values))
		result, covariance = scipy.optimize.curve_fit(linear, x_values, y_values, guess, y_errors)
		print "Linear result: ", result
		errors = numpy.sqrt(numpy.diag(covariance))
		print "Errors:", errors
		m = result[0]
		c = result[1]
		me = errors[0]
		ce = errors[1]
		
		yFit = [linear(x, m, c) for x in spectrum.wavelengths]
		
		ppgplot.pgsci(3)
		ppgplot.pgsls(2)
		ppgplot.pgline(spectrum.wavelengths, yFit)
		ppgplot.pgsls(1)
		

		# d7165
		wavelength, flux = spectrum.getSubsetByWavelength(7140, 7190)
		deltawavelengths = []
		d7165 = 0
		norm = 0
		for i in range(len(wavelength)-1):
			delta = wavelength[i+1] - wavelength[i]
			deltawavelengths.append(delta)
		deltawavelengths.append(delta)
		
		for w, f, dl in zip(wavelength, flux, deltawavelengths):
			d7165+= (linear(w, m, c) - f) * dl/w
			norm+= dl/w
			# print w, f, linear(w, m, c), dl, d7165, norm
			ppgplot.pgsci(4)
			ppgplot.pgline([w, w], [linear(w, m, c), f])
		d7165 = d7165/norm
		print "d7165:", d7165

		# d7665
		wavelength, flux = spectrum.getSubsetByWavelength(7640, 7690)
		deltawavelengths = []
		d7665 = 0
		norm = 0
		for i in range(len(wavelength)-1):
			delta = wavelength[i+1] - wavelength[i]
			deltawavelengths.append(delta)
		deltawavelengths.append(delta)
		
		for w, f, dl in zip(wavelength, flux, deltawavelengths):
			d7665+= (linear(w, m, c) - f) * dl/w
			norm+= dl/w
			print w, f, linear(w, m, c), dl, d7665, norm
			ppgplot.pgsci(4)
			ppgplot.pgline([w, w], [linear(w, m, c), f])
		d7665 = d7665/norm
		print "d7665:", d7665

		# jNa
		wavelength, flux = spectrum.getSubsetByWavelength(8170, 8220)
		deltawavelengths = []
		jNa = 0
		norm = 0
		for i in range(len(wavelength)-1):
			delta = wavelength[i+1] - wavelength[i]
			deltawavelengths.append(delta)
			print delta
		deltawavelengths.append(delta)
		
		for w, f, dl in zip(wavelength, flux, deltawavelengths):
			jNa+= (linear(w, m, c) - f) * 3E8*dl/(w*w)
			# print w, f, linear(w, m, c), dl, d7665, norm
			ppgplot.pgsci(4)
			ppgplot.pgline([w, w], [linear(w, m, c), f])
		print "j Na:", jNa


		print "D ratio:", d7665/d7165
		print "TiO 7165:", d7165 / linear(7500, m, c)
		print "TiO 7665:", d7665 / linear(7500, m, c)
		print "TiO ratio:", d7665 / d7165
	sys.exit()
