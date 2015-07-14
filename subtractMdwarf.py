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

	parser = argparse.ArgumentParser(description='Loads a series of spectra from the JSON files. Subtracts another spectrum from each one.')
	parser.add_argument('inputFiles', type=str, nargs='+', help='JSON files containing the spectra')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	parser.add_argument('--list', action='store_true', help='Specify this option if the input file is actually a list of input files.')
	parser.add_argument('--model', type=str, help='Optional spectrum (to subtract). Also a JSON file.')
	 
	arg = parser.parse_args()
	print arg
	
	
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
		modelSpectrum.objectName = "M-dwarf spectrum"
		modelSpectrum.setData(data['wavelength'], data['flux'])
		tioSpectrum = copy.deepcopy(modelSpectrum)
		tioSpectrum.trimWavelengthRange(7000, 7500)
		tioArea = tioSpectrum.integrate()
		print "tioArea", tioArea

		
	mainPGPlotWindow = ppgplot.pgopen('/xs')	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	yUpper = 2.5
	yLower = -1.0
	
	newSpectra = []
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
	
		subSpectrum = copy.deepcopy(spectrum)
		newLength = subSpectrum.trimWavelengthRange(7000, 7500)
		area = subSpectrum.integrate()
		print "area 7000 - 7500", area
		referenceSpectrum = copy.deepcopy(modelSpectrum)
		# referenceSpectrum.divide(tioArea)
		# referenceSpectrum.divide(1/area)
		checkReference = copy.deepcopy(referenceSpectrum)
		newLength = checkReference.trimWavelengthRange(7000, 7500)
		area = checkReference.integrate()
		print "Reference area:", area
		
		
		ppgplot.pgsci(2)
		ppgplot.pgline(referenceSpectrum.wavelengths, referenceSpectrum.flux)
		
		subtractedSpectrum = copy.deepcopy(spectrum)
		subtractedSpectrum.subtractSpectrum(referenceSpectrum)
		ppgplot.pgsci(3)
		ppgplot.pgline(subtractedSpectrum.wavelengths, subtractedSpectrum.flux)
		
		newSpectra.append(subtractedSpectrum)
		
		
	yUpper = 2.5
	yLower = -1.0
	ppgplot.pgask(True)
	for spectrum in newSpectra:
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
		ppgplot.pgsls(2)
		ppgplot.pgline([lowerWavelength, upperWavelength], [0, 0])
		ppgplot.pgsls(1)
		
		filename = "subtracted_%f.json"%ephemeris.getPhase(spectrum.HJD)
		spectrum.writeToJSON(filename)
		
	
