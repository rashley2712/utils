#!/usr/bin/env python
import sys
import numpy, math
import argparse
import loadingSavingUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot
import generalUtils

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a list of spectra and plots a trail with PGPLOT.')
	parser.add_argument('inputFiles', type=str, nargs='+', help='JSON files containing the spectra')
	parser.add_argument('-e', type=str, help='Optional ephemeris file')
	parser.add_argument('--list', action='store_true', help='Specify this option if the input file is actually a list of input files.')
	parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
	# parser.add_argument('--stacked', action='store_true', help='Specify this option to perform a stacked plot.')
	parser.add_argument('--title', type=str, help='Title for the plot. Otherwise title will be generated from data in the .JSON file.')
	parser.add_argument('--lower', type=float, help='[optional] lower wavelength of the plot.')
	parser.add_argument('--upper', type=float, help='[optional] upper wavelength of the plot.')
	 
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
	
	# Create bitmap for trail plot
	ySize = numSpectra
	xSizeArray = [len(s.wavelengths) for s in spectra]
	xSize = max(xSizeArray)
	print "Bitmap for trails size: (%d, %d)"%(xSize, ySize)
	trailArray = []
	for index, s in enumerate(spectra):
		print index, len(numpy.copy(s.flux))
		trailArray.append(numpy.array(s.flux))
		print numpy.shape(trailArray)
	
	trailBitmap = numpy.copy(trailArray)
	print numpy.shape(trailBitmap)

	mainPGPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(True)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	pgPlotTransform = [0, 1, 0, 0, 0, .1]
	lowerWavelength = min(spectrum.wavelengths)
	upperWavelength = max(spectrum.wavelengths)
	ppgplot.pgenv(0, xSize, 0, ySize, 0, 0)
	ppgplot.pggray(generalUtils.percentiles(trailBitmap, 20, 99), 0, xSize-1 , 0, ySize-1 , 255, 0, pgPlotTransform)
	
	ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "phase [%s]")
		
		
	
