#!/usr/bin/env python
import sys
import numpy, math
import argparse
import loadingSavingUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a series of spectra from the JSON files and joins them in the wavelength dimension.')
	parser.add_argument('files1', type=str, help='List of JSON files containing the first set of spectra')
	parser.add_argument('files2', type=str, help='List of JSON files containing the second set of spectra')
	parser.add_argument('--device', type=str, help='[Optional] PGPLOT device. Defaults to "/xs".')
	parser.add_argument('-e', type=str, help='Ephemeris data file')
	 
	arg = parser.parse_args()
	print arg
	if arg.device!=None:
		device = arg.device
	else:
		device = "/xs"
	
	if arg.e!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.e)
		print ephemeris
	else:
		hasEphemeris = False

	spectra1 = []
	
	filenames = []
	# Load the first set of spectra.
	filename = arg.files1
	fileList = open(filename, 'r')
	for line in fileList:
		filenames.append(str(line.strip()))
	
	print "filenames are:", filenames
	
	spectra1 = []
	for fileIndex, f in enumerate(filenames):
		spectrum = spectrumClasses.spectrumObject()
		spectrum.loadFromJSON(f)
		print "Loaded %s, contains %s."%(f, spectrum.objectName)
		spectra1.append(spectrum)
		
	numSpectra1 = len(spectra1)
	
	filenames = []
	# Load the first set of spectra.
	filename = arg.files2
	fileList = open(filename, 'r')
	for line in fileList:
		filenames.append(str(line.strip()))
	
	spectra2 = []
	for fileIndex, f in enumerate(filenames):
		spectrum = spectrumClasses.spectrumObject()
		spectrum.loadFromJSON(f)
		print "Loaded %s, contains %s."%(f, spectrum.objectName)
		spectra2.append(spectrum)
		
	numSpectra2 = len(spectra2)

		
	mainPGPlotWindow = ppgplot.pgopen(device)	
	ppgplot.pgask(True)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	
	if numSpectra1 != numSpectra2:
		print "We need an equal number of spectra in each list."
		print numSpectra1, numSpectra2
		sys.exit()
		
	for s1, s2 in zip(spectra1, spectra2):
		combinedSpectrum = copy.deepcopy(s1)
		combinedSpectrum.appendNewData(s2)
	
		ppgplot.pgsci(1)
		lowerWavelength = min(combinedSpectrum.wavelengths)
		upperWavelength = max(combinedSpectrum.wavelengths)
		lowerFlux = min(combinedSpectrum.flux)
		upperFlux = max(combinedSpectrum.flux)
		lowerFlux = 0
		upperFlux = 2.5
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
		
		ppgplot.pgsci(2)
		ppgplot.pgline(s1.wavelengths, s1.flux)
		ppgplot.pgsci(3)
		ppgplot.pgline(s2.wavelengths, s2.flux)
		
		
		ppgplot.pgsci(1)
		
		if hasEphemeris:
			ppgplot.pglab("wavelength", "flux", "%s [%f]"%(s1.objectName, ephemeris.getPhase(s1.HJD)))
		else:
			ppgplot.pglab("wavelength", "flux", s1.objectName)
	
		if hasEphemeris:
			filename = "combined_%f.json"%(ephemeris.getPhase(combinedSpectrum.HJD))
			combinedSpectrum.writeToJSON(filename)
