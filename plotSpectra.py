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

def quad(x, a1, a2, a3):
	y = a1 * x * x + a2 *x + a3
	return y

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a series of spectra from the JSON files. Plots them using Matplotlib.')
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

	if arg.model != None:	
		columnNames, data = loadingSavingUtils.loadNewCSV(arg.model)

		modelSpectrum = spectrumClasses.spectrumObject()
		modelSpectrum.objectName = "M-dwarf model"
		modelSpectrum.setData(data['wavelength'], data['flux'])


	for s in spectra:
		if hasEphemeris:
			uniqueString = "%1.3f"%ephemeris.getPhase(s.HJD)
		else:
			uniqueString = "%f"%s.HJD
		filename = s.objectName + "_" + uniqueString + ".json"
		s.writeToJSON(filename)
	
	test = spectrumClasses.spectrumObject()
	test.loadFromJSON("OT0711_0.815.json")

	print test.ra


	matplotlib.pyplot.figure(figsize=(20, 4*numSpectra + 1))
	offset = 1.0
	for index, spectrum in enumerate(spectra):
		if arg.model!=None:
			spectrum.subtractSpectrum(modelSpectrum)
		yValues = [y + offset * index for y in spectrum.getFlux()]
		matplotlib.pyplot.plot(spectrum.getWavelengths(), yValues, drawstyle = 'steps', color = 'k')
		if hasEphemeris:
			phase = ephemeris.getPhase(spectrum.HJD)
			labelX = 5550
			labelY = spectrum.getNearestFlux(labelX) + offset*index + offset/5.0
			matplotlib.pyplot.text(labelX, labelY, 'phase: %1.2f'%(phase), fontsize=20)
	
	matplotlib.pyplot.yticks(fontsize = 22)
	matplotlib.pyplot.ylabel(r"F$_{\nu}$ (mJy) with offset of %1.1f applied."%offset, size = 22)
	matplotlib.pyplot.xticks(fontsize = 22)
	matplotlib.pyplot.xlabel(r"Wavelength $(\AA)$", size = 22)

	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show(block = False)
	fig.savefig('spectra.eps',dpi=100, format='eps')
	fig.savefig('spectra.png',dpi=200, format='png')

	matplotlib.pyplot.figure(figsize=(16, 4*numSpectra + 1))
	offset = 0.7
	for index, spectrum in enumerate(spectra):
		(wavelengths, flux) = spectrum.getSubsetByWavelength(8000, 8400)
		yValues = [3*y + offset * index for y in flux]
		matplotlib.pyplot.plot(wavelengths, yValues, color = 'k', drawstyle='steps')
		if hasEphemeris:
			phase = ephemeris.getPhase(spectrum.HJD)
			labelX = 8010
			labelY = spectrum.getNearestFlux(labelX) + offset*index + offset/5.0
			matplotlib.pyplot.text(labelX, labelY, 'phase: %1.2f'%(phase), fontsize=15)
	
	matplotlib.pyplot.ylabel(r"F$_{\nu}$ (mJy)", size = 16)
	matplotlib.pyplot.xlabel(r"Wavelength $(\AA)$", size = 16)

	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show(block= True)
	fig.savefig('sodium8190.eps',dpi=100, format='eps')
	fig.savefig('sodium8190.png',dpi=200, format='png')

	spectrum = copy.deepcopy(spectra[8])
	numSpectra = 1
	matplotlib.pyplot.figure(figsize=(20, 4*numSpectra + 1))
	offset = 0
	yValues = [y + offset * index for y in spectrum.getFlux()]
	matplotlib.pyplot.plot(spectrum.getWavelengths(), yValues, drawstyle = 'steps', color = 'k')
	if hasEphemeris:
		phase = ephemeris.getPhase(spectrum.HJD)
		labelX = 5550
		labelY = spectrum.getNearestFlux(labelX) + offset*index + 1./5.0
		matplotlib.pyplot.text(labelX, labelY, 'phase: %1.2f'%(phase), fontsize=15)
	
	matplotlib.pyplot.ylabel(r"F$_{\nu}$ (mJy)", size = 16)
	matplotlib.pyplot.xlabel(r"Wavelength $(\AA)$", size = 16)

	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show(block = False)
	fig.savefig('spectrum.eps',dpi=100, format='eps')
	fig.savefig('spectrum.png',dpi=200, format='png')

	# Now snip out part of the spectrum
	spectrum = copy.deepcopy(spectra[10])
	spectrum.trimWavelengthRange(6000, 6850)
	spectrum.snipWavelengthRange(6540, 6600)
	
	fluxValues = spectrum.getFlux()
	wavelengthValues = spectrum.getWavelengths()
	
	x_values = wavelengthValues
	y_values = fluxValues
	y_errors = numpy.ones(len(x_values))
	# print x_values, y_values, y_errors
	a1 = 0.0
	a2 = 0.0
	a3 = 0.0 
	guess = numpy.array([a1, a2, a3])
	print "initial guess", guess
	result, covariance = scipy.optimize.curve_fit(quad, x_values, y_values, guess, y_errors)
	print "Quadratic result: ", result
	errors = numpy.sqrt(numpy.diag(covariance))
	print "Errors:", errors
	a1 = result[0]
	a2 = result[1]
	a3 = result[2]
	a1e = errors[0]
	a2e = errors[1]
	a3e = errors[2]
	
	y_model = [quad(x, a1, a2, a3) for x in x_values]
	peak = -1.0*a2 / (2.0 * a1)	
	peakError = peak * math.sqrt( (a1e/a1)**2 + (a2e/a2)**2)
	print "quadratic peaks at:", peak, peakError
	
	numSpectra = 1
	matplotlib.pyplot.figure(figsize=(20, 4*numSpectra + 1))
	matplotlib.pyplot.plot(wavelengthValues, fluxValues, drawstyle = 'steps', color = 'k')
	matplotlib.pyplot.plot(wavelengthValues, y_model, color='b', linestyle='dashed')
	matplotlib.pyplot.ylabel(r"F$_{\nu}$ (mJy)", size = 16)
	matplotlib.pyplot.xlabel(r"Wavelength $(\AA)$", size = 16)
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('hump1.eps',dpi=100, format='eps')
	fig.savefig('hump1.png',dpi=200, format='png')

	# Now snip out part of the spectrum
	spectrum = copy.deepcopy(spectra[10])
	spectrum.trimWavelengthRange(7700, 8500)
	spectrum.snipWavelengthRange(8170, 8200)
	
	fluxValues = spectrum.getFlux()
	wavelengthValues = spectrum.getWavelengths()
	
	x_values = wavelengthValues
	y_values = fluxValues
	y_errors = numpy.ones(len(x_values))
	# print x_values, y_values, y_errors
	a1 = 0.0
	a2 = 0.0
	a3 = 0.0 
	guess = numpy.array([a1, a2, a3])
	print "initial guess", guess
	result, covariance = scipy.optimize.curve_fit(quad, x_values, y_values, guess, y_errors)
	print "Quadratic result: ", result
	errors = numpy.sqrt(numpy.diag(covariance))
	print "Errors:", errors
	a1 = result[0]
	a2 = result[1]
	a3 = result[2]
	a1e = errors[0]
	a2e = errors[1]
	a3e = errors[2]
	
	y_model = [quad(x, a1, a2, a3) for x in x_values]
	peak = -1.0*a2 / (2.0 * a1)	
	peakError = peak * math.sqrt( (a1e/a1)**2 + (a2e/a2)**2) 
	print "quadratic peaks at:", peak, peakError
	
	numSpectra = 1
	matplotlib.pyplot.figure(figsize=(20, 4*numSpectra + 1))
	matplotlib.pyplot.plot(wavelengthValues, fluxValues, drawstyle = 'steps', color = 'k')
	matplotlib.pyplot.plot(wavelengthValues, y_model, color='b', linestyle='dashed')
	matplotlib.pyplot.ylabel(r"F$_{\nu}$ (mJy)", size = 16)
	matplotlib.pyplot.xlabel(r"Wavelength $(\AA)$", size = 16)
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show()
	fig.savefig('hump2.eps',dpi=100, format='eps')
	fig.savefig('hump2.png',dpi=200, format='png')


	# write all spectra to CSV
	for spectrum in spectra:
		phase = ephemeris.getPhase(spectrum.HJD)
		filename = "phase%f.csv"%phase
		spectrum.writeCSV(filename)
	
