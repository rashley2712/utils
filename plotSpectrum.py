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

	parser = argparse.ArgumentParser(description='Loads a spectrum JSON file and plots it with PGPLOT.')
	parser.add_argument('inputFiles', type=str, nargs='+', help='JSON files containing the spectra')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	parser.add_argument('--list', action='store_true', help='Specify this option if the input file is actually a list of input files.')
	parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
	parser.add_argument('--stacked', action='store_true', help='Specify this option to perform a stacked plot.')
	parser.add_argument('--title', type=str, help='Title for the plot. Otherwise title will be generated from data in the .JSON file.')
	parser.add_argument('--lower', type=float, help='[optional] lower wavelength of the plot.')
	parser.add_argument('--upper', type=float, help='[optional] upper wavelength of the plot.')
	parser.add_argument('-n', '--normalise', action='store_true', help='Perform a normalise on the spectra. Mean value will be taken from the first spectrum between the ''-nu'' ''-nl'' wavelengths.')
	parser.add_argument('-nu', type=float, help='Upper wavelength of the spectrum for the normalisation average. Required if ''-n'' is specified.')
	parser.add_argument('-nl', type=float, help='Lower wavelength of the spectrum for the normalisation average. Required if ''-n'' is specified.')
	
	 
	arg = parser.parse_args()
	print arg
	
	if arg.normalise:
		if arg.nu is None or arg.nl is None:
			print "We require a '-nu' value to perform the normalise function."
			sys.exit()
	
	if arg.e!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.e)
		print ephemeris
	else:
		hasEphemeris = False

	
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
	
	
	spectra = []
	for fileIndex, f in enumerate(filenames):
		spectrum = spectrumClasses.spectrumObject()
		spectrum.loadFromJSON(f)
		print "Loaded %s, contains %s."%(f, spectrum.objectName)
		if hasEphemeris:
			phase = ephemeris.getPhase(spectrum.HJD)
			spectrum.phase = phase
		spectra.append(spectrum)
		
	numSpectra = len(spectra)
	
	if hasEphemeris:
		# Sort the spectra by their phase
		spectra = sorted(spectra, key=lambda object: object.phase, reverse = False)
			
	numSpectra = len(spectra)
	if numSpectra>1:
		print "%d spectra have been loaded."%numSpectra
		
	if (arg.upper != None) and (arg.lower != None):
		for s in spectra:
			s.trimWavelengthRange(arg.lower, arg.upper)	
	
	if arg.normalise:
		# Perform the normalisation across all spectra
		referenceSpectrum = spectra[0]
		normalConstant = referenceSpectrum.integrate((arg.nl, arg.nu))
		print "Normalisation constant:", normalConstant
	
		for index in range(1, len(spectra)):
			s = spectra[index]
			normalVal = s.integrate((arg.nl, arg.nu))
			print "Normalisation value:", normalVal, normalConstant
			spectra[index].divide(normalVal/normalConstant)
			normalVal = s.integrate((arg.nl, arg.nu))
			print "New Normalisation value:", normalVal, normalConstant
	
	
	
	if not arg.stacked:
		mainPGPlotWindow = ppgplot.pgopen(arg.device)	
		ppgplot.pgask(True)
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		yUpper = 2.5
		yLower = -0.5
		for spectrum in spectra:
			ppgplot.pgsci(1)
			lowerWavelength = min(spectrum.wavelengths)
			upperWavelength = max(spectrum.wavelengths)
			lowerFlux = min(spectrum.flux)
			upperFlux = max(spectrum.flux)
			ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
			ppgplot.pgbin(spectrum.wavelengths, spectrum.flux)
			if hasEphemeris:
				ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "flux [%s]"%spectrum.fluxUnits, "%s [%f]"%(spectrum.objectName, spectrum.phase))
			else:
				ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "flux [%s]"%spectrum.fluxUnits, "%s [%s]"%(spectrum.objectName, spectrum.loadedFromFilename))
		
		
	
	if arg.stacked:
		mainPGPlotWindow = ppgplot.pgopen(arg.device)	
		ppgplot.pgask(True)
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		yLower = -0.5
		offset = 2.0
		yUpper = numSpectra * offset
		lowerWavelength = min(spectrum.wavelengths)
		upperWavelength = max(spectrum.wavelengths)
		ppgplot.pgpap(6.18, 1.618)
		ppgplot.pgenv(lowerWavelength, upperWavelength, yLower, yUpper, 0, 0)
		for index, spectrum in enumerate(spectra):
			ppgplot.pgsci(1)
			flux = [ f + offset*index for f in spectrum.flux]
			ppgplot.pgbin(spectrum.wavelengths, flux)
			plotx = 8600
			ploty = offset*index + spectrum.getNearestFlux(plotx) + offset/5
			if hasEphemeris:
				ppgplot.pgptxt(plotx, ploty, 0, 0, "%1.2f"%(spectrum.phase)) 
			if arg.title!=None:
				ppgplot.pglab("wavelength", "flux", arg.title)
			else:
				ppgplot.pglab("wavelength", "flux", spectrum.objectName)
		
