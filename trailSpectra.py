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
	parser.add_argument('--phasebins', type=int, default=50, help='[optional] number of phase bins to use. Default is 50.')
	parser.add_argument('-n', '--normalise', action='store_true', help='Perform a normalise on the spectra. Mean value will be taken from the first spectrum between the ''-nu'' ''-nl'' wavelengths.')
	parser.add_argument('-nu', type=float, help='Upper wavelength of the spectrum for the normalisation average. Required if ''-n'' is specified.')
	parser.add_argument('-nl', type=float, help='Lower wavelength of the spectrum for the normalisation average. Required if ''-n'' is specified.')
	
	 
	arg = parser.parse_args()
	print arg
	
	if arg.normalise:
		if arg.nu is None or arg.nl is None:
			print "We require a '-nu' and a '-nl' value to perform the normalise function."
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
			print spectrum.HJD, spectrum.phase
		spectra.append(spectrum)
		
	numSpectra = len(spectra)
	
	if hasEphemeris:
		# Sort the spectra by their phase
		spectra = sorted(spectra, key=lambda object: object.phase, reverse = False)
			
	numSpectra = len(spectra)
	if numSpectra>1:
		print "%d spectra have been loaded."%numSpectra
	
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
	
	
	if (arg.upper != None) and (arg.lower != None):
		numPoints = -1
		wavelengthTemplate = []
		for s in spectra:
			s.trimWavelengthRange(arg.lower, arg.upper)	
			if not numPoints==-1:
				s.resample(wavelengthTemplate)
			else:
				numPoints = len(s.wavelengths)
				wavelengthTemplate = s.wavelengths
	
	
	

	# sys.exit()
	
	# Create bitmap for trail plot
	if not hasEphemeris:
		ySize = numSpectra
		xSizeArray = [len(s.wavelengths) for s in spectra]
		xSize = max(xSizeArray)
		print "Bitmap for trails size: (%d, %d)"%(xSize, ySize)
		trailArray = []
		for index, s in enumerate(spectra):
			trailArray.append(numpy.array(s.flux))
	
		trailBitmap = numpy.copy(trailArray)

	else:
		numPhaseBins = arg.phasebins
		ySize = numPhaseBins * 2
		xSizeArray = [len(s.wavelengths) for s in spectra]
		xSize = max(xSizeArray)
		trailArray = []
		for index in range(ySize):
			blankSpectrum = numpy.zeros(xSize)
			trailArray.append(blankSpectrum)
		for s in spectra:
			phase = s.phase
			phaseBin = int(phase*numPhaseBins)
			print phase, phaseBin
			existingSpectrum = trailArray[phaseBin]
			newSpectrum = (existingSpectrum + s.flux )/2
			trailArray[phaseBin] = newSpectrum
			# Also add the spectrum to the phase+1 bin
			phase = s.phase + 1.0
			phaseBin = int(phase*float(numPhaseBins))
			print phase, phaseBin
			existingSpectrum = trailArray[phaseBin]
			newSpectrum = (numpy.array(existingSpectrum) + numpy.array(s.flux)) / 2
			trailArray[phaseBin] = newSpectrum
			


		trailBitmap = numpy.copy(trailArray)
		

	mainPGPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(True)
	ppgplot.pgpap(5.0, 2)
	# pgPlotTransform = [0, 1, 0, 0, 0, 1]
	spectrum = spectra[0] 
	xScale = (max(spectrum.wavelengths) - min(spectrum.wavelengths)) / xSize
	if hasEphemeris:
		yScale = numPhaseBins
	else:	
		yScale = len(spectra)
		
	pgPlotTransform = [min(spectrum.wavelengths), xScale, 0, 0, 0, 1/float(yScale)]
	lowerWavelength = min(spectrum.wavelengths)
	upperWavelength = max(spectrum.wavelengths)
	ppgplot.pgenv( min(spectrum.wavelengths), max(spectrum.wavelengths), 0, 2, 0, 0)
	ppgplot.pggray(generalUtils.percentiles(trailBitmap, 20, 99), 0, xSize-1 , 0, ySize-1 , 255, 0, pgPlotTransform)
	# ppgplot.pggray(trailBitmap, 0, xSize-1 , 0, ySize-1 , numpy.max(trailBitmap), numpy.min(trailBitmap), pgPlotTransform)
	
	epochs = [s.HJD for s in spectra]
	startHJD = min(epochs)
	endHJD = max(epochs)
	if arg.title is None:
		title = ""
	else: 
		title = arg.title
	if hasEphemeris: 
		# ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "Phase", "%s - %s"%(str(startHJD), str(endHJD)))
		ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "Phase", title)
	else:
		#ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "Spectrum number", "%s - %s"%(str(spectra[0].HJD), str(spectra[-1].HJD)))
		ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "Spectrum number", title)
		
		
		
	
