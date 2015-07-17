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

	parser = argparse.ArgumentParser(description='Loads a series of spectra from the JSON files and tries to fit cyclotron humps.')
	parser.add_argument('files', type=str, help='List of JSON files containing the first set of spectra')
	parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
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

	spectra = []
	
	filenames = []
	# Load the first set of spectra.
	filename = arg.files
	fileList = open(filename, 'r')
	for line in fileList:
		filenames.append(str(line.strip()))
	
	for fileIndex, f in enumerate(filenames):
		spectrum = spectrumClasses.spectrumObject()
		spectrum.loadFromJSON(f)
		print "Loaded %s, contains %s."%(f, spectrum.objectName)
		spectra.append(spectrum)
		
	
	mainPGPlotWindow = ppgplot.pgopen(device)	
	ppgplot.pgask(True)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	
		
	for spectrum in spectra:
		
		ppgplot.pgsci(1)
		lowerWavelength = min(spectrum.wavelengths)
		upperWavelength = max(spectrum.wavelengths)
		lowerFlux = min(spectrum.flux)
		upperFlux = max(spectrum.flux)
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
		
		ppgplot.pgline(spectrum.wavelengths, spectrum.flux)
		
		ppgplot.pgsci(1)
		
		if hasEphemeris:
			ppgplot.pglab("wavelength", "flux", "%s [%f]"%(spectrum.objectName, ephemeris.getPhase(spectrum.HJD)))
		else:
			ppgplot.pglab("wavelength", "flux", spectrum.objectName)
	
	
