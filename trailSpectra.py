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
			print spectrum.HJD, spectrum.phase
		spectra.append(spectrum)
		
	numSpectra = len(spectra)
	
	if hasEphemeris:
		# Sort the spectra by their phase
		spectra = sorted(spectra, key=lambda object: object.phase, reverse = False)
			
	numSpectra = len(spectra)
	if numSpectra>1:
		print "%d spectra have been loaded."%numSpectra
	
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
		numPhaseBins = 40
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
			newSpectrum = (existingSpectrum + s.flux )/2
			trailArray[phaseBin] = newSpectrum
			


		trailBitmap = numpy.copy(trailArray)
		

	mainPGPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(True)
	# pgPlotTransform = [0, 1, 0, 0, 0, 1]
	spectrum = spectra[0] 
	xScale = (max(spectrum.wavelengths) - min(spectrum.wavelengths)) / xSize
	pgPlotTransform = [min(spectrum.wavelengths), xScale, 0, 0, 0, 1/float(numPhaseBins)]
	lowerWavelength = min(spectrum.wavelengths)
	upperWavelength = max(spectrum.wavelengths)
	ppgplot.pgenv( min(spectrum.wavelengths), max(spectrum.wavelengths), 0, 2, 0, 0)
	ppgplot.pggray(generalUtils.percentiles(trailBitmap, 20, 99), 0, xSize-1 , 0, ySize-1 , 255, 0, pgPlotTransform)
	
	epochs = [s.HJD for s in spectra]
	startHJD = min(epochs)
	endHJD = max(epochs)
	if hasEphemeris: 
		ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "Phase", "%s - %s"%(str(startHJD), str(endHJD)))
	else:
		ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "Spectrum number", "%s - %s"%(str(spectra[0].HJD), str(spectra[-1].HJD)))
		
		
		
	
