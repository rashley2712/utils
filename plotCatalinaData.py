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
import trm.dnl.molly
import spectrumClasses, timeClasses
import scipy.optimize
import copy

def quad(x, a1, a2, a3):
	y = a1 * x * x + a2 *x + a3
	return y

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a series of spectra that were saved from Molly. Plots them using Matplotlib.')
	parser.add_argument('inputFiles', type=str, nargs='+', help='Molly files containing the spectra')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	 
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
	
	for fileIndex, f in enumerate(arg.inputFiles):
	
		mollyFile = trm.dnl.molly.rmolly(f)
	
		
		for index, r in enumerate(mollyFile):
			data = r.toHDU()
			wavelengths = r.x.data
			flux = r.y.data
			head = r.head
			if r.has_mask: print "Mask found"
			spectrum = spectrumClasses.spectrumObject()
			npoints = spectrum.setData(wavelengths, flux)
			targetName = spectrum.parseHeaderInfo(head)
			print "Parsed headers of", targetName
			print r.oneLine()
			if fileIndex == 0:
				spectra.append(spectrum)
			else:
				#spectra[index].appendDataAtNewWavelengths(wavelengths, flux)
				spectra[index].appendNewData(spectrum)
		
	
		
	numSpectra = len(spectra)
	matplotlib.pyplot.figure(figsize=(20, 4*numSpectra + 1))
	offset = 1.2
	for index, spectrum in enumerate(spectra):
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

	numSpectra = len(spectra)
	matplotlib.pyplot.figure(figsize=(16, 4*numSpectra + 1))
	offset = 1
	for index, spectrum in enumerate(spectra):
		(wavelengths, flux) = spectrum.getSubsetByWavelength(6400, 6800)
		yValues = [y + offset * index for y in flux]
		matplotlib.pyplot.plot(wavelengths, yValues, drawstyle = 'steps', color = 'k')
		if hasEphemeris:
			phase = ephemeris.getPhase(spectrum.HJD)
			labelX = 6450
			labelY = spectrum.getNearestFlux(labelX) + offset*index + offset/5.0
			matplotlib.pyplot.text(labelX, labelY, 'phase: %1.2f'%(phase), fontsize=15)
	
	matplotlib.pyplot.ylabel(r"F$_{\nu}$ (mJy)", size = 16)
	matplotlib.pyplot.xlabel(r"Wavelength $(\AA)$", size = 16)

	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.show(block= False)
	fig.savefig('Halpha_zoom.eps',dpi=100, format='eps')
	fig.savefig('Halpha_zoom.png',dpi=200, format='png')


	spectrum = copy.deepcopy(spectra[10])
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
