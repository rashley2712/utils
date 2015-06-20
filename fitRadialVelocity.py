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
import ppgplot


def poly(x, n, a):
	y = 0;
	coeffs = reversed(a)
	for index in range(n):
		y+= (x**index) * coeffs[index] 
	return y

def quad(x, a1, a2, a3):
	y = a1 * x * x + a2 *x + a3
	return y

def gaussian(x, a0, a1, a2, a3):
	y = a0 + a1 * numpy.exp(-.5 * ((x-a2)/a3)**2)
	return y

def doubleGaussian(x, a0, a1, a2):
	s = 11.0
	w = 3.4
	y = a0 + a1 * numpy.exp(-.5 * ((x-a2)/w)**2) + a1 * numpy.exp(-.5 * (((x-(a2+s))/w)**2) )
	print a0, a1, a2
	return y

def sinewave(x, a0, a1, a2):
	omega = 2*math.pi
	y = a0 + a1 * numpy.sin( omega * (x + a2))
	return y

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a series of spectra that were saved from Molly. Tries to fit radial velocity profile to the spectra.')
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
				spectra[index].appendNewData(spectrum)
				
	numSpectra = len(spectra)

	print "%d spectra have been loaded."%numSpectra

	mainPGPlotWindow = ppgplot.pgopen('/xs')	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgask(False)
	ppgplot.pglab("wavelength", "flux", "spectrum")
	phases = []
	rvs = []
	wvshifts = []
	wvshiftErrors = []
	previousWavelengthFit = 8173.0
	for spectrum in spectra:
		spectrum.trimWavelengthRange(8000, 8400)
		lowerWavelength = min(spectrum.wavelengths)
		upperWavelength = max(spectrum.wavelengths)
		lowerFlux = min(spectrum.flux)
		upperFlux = max(spectrum.flux)
		ppgplot.pgsci(1)
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
		ppgplot.pgline(spectrum.wavelengths, spectrum.flux)
		if hasEphemeris:
			phase = ephemeris.getPhase(spectrum.HJD)
			labelX = 8100
			labelY = spectrum.getNearestFlux(labelX) + 0.5
			print labelX, labelY
			ppgplot.pglab("Wavelength", "Flux", "Phase: %f"%phase)

		# Now remove the 8190 doublet from the data set
		continuum = copy.deepcopy(spectrum)
		continuum.snipWavelengthRange(8130, 8250)
		ppgplot.pgsci(3)
		ppgplot.pgline(continuum.wavelengths, continuum.flux)
		
		x_values = continuum.wavelengths
		y_values = continuum.flux
		y_errors = numpy.ones(len(x_values))
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
		
		y_fit = [quad(x, a1, a2, a3) for x in spectrum.wavelengths]
		ppgplot.pgsci(2)
		ppgplot.pgline(spectrum.wavelengths, y_fit)
		
		flatSpectrumFlux = [data/fit for (data, fit) in zip(spectrum.flux, y_fit)]
		flatSpectrum = copy.deepcopy(spectrum)
		flatSpectrum.setData(spectrum.wavelengths, flatSpectrumFlux)

		ppgplot.pgsci(5)
		ppgplot.pgenv(lowerWavelength, upperWavelength, min(flatSpectrum.flux), max(flatSpectrum.flux), 0, 0)
		ppgplot.pgline(flatSpectrum.wavelengths, flatSpectrum.flux)
		ppgplot.pglab("Wavelength", "Flux", "Name: %s Phase: %f"%(flatSpectrum.objectName, phase))

		flatSpectrum.writeCSV("Na_i_%s_%f.csv"%(flatSpectrum.objectName, phase))
		subSpectrum = copy.deepcopy(flatSpectrum)
		subSpectrum.trimWavelengthRange(8100, 8290)


		a0 = 1.0
		a1 = -0.3
		a2 = previousWavelengthFit
		gaussianX = numpy.arange(min(subSpectrum.wavelengths), max(subSpectrum.wavelengths), 1)
		gaussianY = doubleGaussian(gaussianX, a0, a1, a2)
		# ppgplot.pgsci(2)
		# ppgplot.pgline(gaussianX, gaussianY)
		guess = numpy.array([a0, a1, a2])
		x_values = subSpectrum.wavelengths
		y_values = subSpectrum.flux
		y_errors = numpy.ones(len(x_values))
		result, covariance = scipy.optimize.curve_fit(doubleGaussian, x_values, y_values, guess, y_errors)
		errors = numpy.sqrt(numpy.diag(covariance))
		print "Result:", result
		print "Errors:", errors
		a0 = result[0]
		a1 = result[1]
		a2 = result[2]
		gaussianX = numpy.arange(min(subSpectrum.wavelengths), max(subSpectrum.wavelengths), 0.1)
		gaussianY = doubleGaussian(gaussianX, a0, a1, a2)
		ppgplot.pgsci(1)
		ppgplot.pgline(gaussianX, gaussianY)
		print "Phase:%f Wavelength:%f [%f]"%(phase, a2, errors[2])
		previousWavelengthFit = a2
		phases.append(phase)
		wvshifts.append(a2)
		wvshiftErrors.append(errors[2])

	# Convert the wavelengths into velocities. 
	velocity = [(w- 8183)/8183 * 3E5 for w in wvshifts]
	velocityError = [ 3E5 /8183 * d for d in wvshiftErrors]
	print "Velocities:", velocity
	print "Velocity errors:", velocityError
	
	# Now fit a sine wave to our plot
	a0 = 0.0  
	a1 = 250.0
	a2 = 0.0
	guess = numpy.array([a0, a1, a2])
	x_values = phases
	y_values = velocity
	y_errors = velocityError
	results, covariance = scipy.optimize.curve_fit(sinewave, x_values, y_values, guess, y_errors)
	print "Sine result:", results
	print "Sine errors:", covariance
	a0 = results[0]
	a1 = results[1]
	a2 = results[2]
	xFit = numpy.arange(0, 1, 0.01)
	yFit = sinewave(xFit, a0, a1, a2)

	mainPGPlotWindow = ppgplot.pgopen('/xs')	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgsci(1)
	ppgplot.pgenv(0, 1, min(velocity)*1.1, max(velocity)*1.1, 0, 0)
	ppgplot.pgpt(phases, velocity)
	ppgplot.pgerrb(2, phases, velocity, velocityError, 0)
	ppgplot.pgerrb(4, phases, velocity, velocityError, 0)
	ppgplot.pgsci(5)
	ppgplot.pgsls(2)
	ppgplot.pgline( [0, 1], [a0, a0])
	ppgplot.pgsls(1)
	ppgplot.pgline(xFit, yFit)
	ppgplot.pgsci(1)
	ppgplot.pglab("phase [phi]", "velocity [km/s]", "%.0f + %.0f*sin(2*pi*[phi + %.4f])"%(a0, a1, a2))
	
	mainPGPlotWindow = ppgplot.pgopen('pgplot.eps/ps')	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgsci(1)
	ppgplot.pgenv(0, 1, min(velocity)*1.1, max(velocity)*1.1, 0, 0)
	ppgplot.pgpt(phases, velocity)
	ppgplot.pgerrb(2, phases, velocity, velocityError, 0)
	ppgplot.pgerrb(4, phases, velocity, velocityError, 0)
	ppgplot.pgsci(5)
	ppgplot.pgsls(2)
	ppgplot.pgline( [0, 1], [a0, a0])
	ppgplot.pgsls(1)
	ppgplot.pgline(xFit, yFit)
	ppgplot.pgsci(1)
	ppgplot.pglab("phase [phi]", "velocity [km/s]", "")
	
	for p, w, e in zip(phases, wvshifts, wvshiftErrors):
		print p, w, e
