#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import loadingSavingUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot

def quad(x, a0, a1, a2):
	y = a2 * x * x + a1 *x + a0
	return y

def gaussian(x, a0, a1, a2, a3):
	y = a0 + a1 * numpy.exp(-.5 * ((x-a2)/a3)**2)
	return y	
	
# Lab wavelengths for Sodium doublet  8183 and 8195
# 8183.2556 and 8194.7905   separation: 11.5349
def doubleGaussian(x, a0, a1, a2):
	global width
	s = 11.5349
	w = 3.4
	w = width
	y = a0 + a1 * numpy.exp(-.5 * ((x-a2)/w)**2) + a1 * numpy.exp(-.5 * (((x-(a2+s))/w)**2) )
	return y

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a spectrum JSON file and fits a double gaussian to the 8190 Na doublet.')
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
	parser.add_argument('--width', type=float, help='FWHM parameter for the Gaussian fit. If specified, this value will be fixed. ')
	parser.add_argument('-o', '--objectname', type=str, help='[Optional] Object name for the output log file.')
	
	
	 
	arg = parser.parse_args()
	# print arg
	
	if arg.width is not None:
		fixWidth = True
		width = arg.width
		print "Fixing line width to %f AA."%width
	else:
		fixWidth = False
		
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

	trimLower = 8150
	trimUpper = 8240
	print "Trimming out the region around 8190AA. [%d, %d]"%(trimLower, trimUpper)
	for s in spectra:
		s.trimWavelengthRange(trimLower, trimUpper)	

	
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
	
	
	
	mainPGPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(True)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	yUpper = 2.5
	yLower = -0.5
	
	fitPGPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(True)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
		
	
	for spectrum in spectra:
		ppgplot.pgslct(mainPGPlotWindow)
		
		ppgplot.pgsci(1)
		lowerWavelength = min(spectrum.wavelengths)
		upperWavelength = max(spectrum.wavelengths)
		lowerFlux = min(spectrum.flux)
		upperFlux = max(spectrum.flux)
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
		ppgplot.pgbin(spectrum.wavelengths, spectrum.flux)
		ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "flux [%s]"%spectrum.fluxUnits, "%s [%s]"%(spectrum.objectName, spectrum.loadedFromFilename))
		
		# Grab the continuum from either side of the spectrum
		
		lowerCut = 8170
		upperCut = 8220
		continuumSpectrum = copy.deepcopy(spectrum)
		continuumSpectrum.snipWavelengthRange(lowerCut, upperCut)
		ppgplot.pgsci(2)
		ppgplot.pgbin(continuumSpectrum.wavelengths, continuumSpectrum.flux)
		
		# Now fit a polynomial to continuum around the doublet 
		a0 = 0.0    	# Constant term  
		a0 = numpy.mean(continuumSpectrum.flux)
		a1 = 0.0		# Linear term
		a2 = 0.0		# Quadratic term
		guess = numpy.array([a0, a1, a2])
		x_values = continuumSpectrum.wavelengths
		y_values = continuumSpectrum.flux
		y_errors = numpy.ones(len(continuumSpectrum.flux))
		results, covariance = scipy.optimize.curve_fit(quad, x_values, y_values, guess, )
		errors = numpy.sqrt(numpy.diag(covariance))
		print "quadratic result:", results
		print "quadratic errors:", errors
		a0 = results[0]
		a1 = results[1]
		a2 = results[2]
		xFit = spectrum.wavelengths
		yFit = quad(numpy.array(xFit), a0, a1, a2)

		ppgplot.pgsci(3)
		ppgplot.pgline(xFit, yFit)
	
		normalisedSpectrum = copy.deepcopy(spectrum)
		for index, w in enumerate(spectrum.wavelengths):			
			# print w, spectrum.flux[index], yFit[index], spectrum.flux[index]/yFit[index]
			normalisedSpectrum.flux[index] = spectrum.flux[index]/yFit[index]
		
		lowerWavelength = min(normalisedSpectrum.wavelengths)
		upperWavelength = max(normalisedSpectrum.wavelengths)
		lowerFlux = min(normalisedSpectrum.flux)
		upperFlux = max(normalisedSpectrum.flux)
	
		ppgplot.pgslct(fitPGPlotWindow)
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
		ppgplot.pgsci(1)
		ppgplot.pgbin(normalisedSpectrum.wavelengths, normalisedSpectrum.flux)
		
		featureSpectrum = copy.deepcopy(normalisedSpectrum)
		featureSpectrum.trimWavelengthRange(lowerCut, upperCut)
		# Now fit a single gaussian to the Na doublet blue line 
		a0 = 1.0    	# Constant term  
		a1 = -0.5		# 'Depth' of the line
		a2 = 8183.0		# Wavelength of the centre of the line
		a3 = 3.0		# Width of the line
		guess = numpy.array([a0, a1, a2, a3])
		x_values = featureSpectrum.wavelengths
		y_values = featureSpectrum.flux
		y_errors = numpy.ones(len(featureSpectrum.flux))
		results, covariance = scipy.optimize.curve_fit(gaussian, x_values, y_values, guess, y_errors)
		errors = numpy.sqrt(numpy.diag(covariance))
		print "gaussian result:", results
		print "gaussian errors:", errors
		a0 = results[0]
		a1 = results[1]
		a2 = results[2]
		a3 = results[3]
		print "Centroid wavelength %f [%f]"%(a2, errors[2])
		print "Width %f [%f]"%(a3, errors[3])
		xFit = spectrum.wavelengths
		yFit = gaussian(numpy.array(xFit), a0, a1, a2, a3)
		
		if not fixWidth:
			width = a3
		centroidWavelength = a2
		depth = a1
		constant = a0

		ppgplot.pgsci(3)
		ppgplot.pgline(xFit, yFit)

		# Now fit the double gaussian with a fixed width and separation
		a0 = constant    			# Constant term  
		a1 = depth					# 'Depth' of the line
		a2 = centroidWavelength		# Wavelength of the centre of the blueward line
		guess = numpy.array([a0, a1, a2])
		x_values = featureSpectrum.wavelengths
		y_values = featureSpectrum.flux
		y_errors = numpy.ones(len(featureSpectrum.flux))
		results, covariance = scipy.optimize.curve_fit(doubleGaussian, x_values, y_values, guess, y_errors)
		errors = numpy.sqrt(numpy.diag(covariance))
		print "double gaussian result:", results
		print "double gaussian errors:", errors
		a0 = results[0]
		a1 = results[1]
		a2 = results[2]
		wavelength = a2
		wavelengthError = errors[2]
		print "Centroid blueward wavelength %f [%f]"%(wavelength, wavelengthError)
		velocity = (wavelength - 8183)/8183 * 3E5 
		velocityError = 3E5 /8183 * wavelengthError
		print "Velocity %f [%f]"%(velocity, velocityError)

		xFit = spectrum.wavelengths
		yFit = doubleGaussian(numpy.array(xFit), a0, a1, a2)
		ppgplot.pgsci(4)
		ppgplot.pgline(xFit, yFit)

		ppgplot.pgsci(5)
		ppgplot.pgsls(2)
		ppgplot.pgline([8183, 8183], [lowerFlux, upperFlux])
		ppgplot.pgsls(1)
		
		if not query_yes_no("Are you happy with the fit?"):
			print "Not saving the value. Fit not good enough"
			sys.exit()
		
		# Now record the RV in a log file
		logLines = []
		oldLines = []
		if arg.objectname is not None:
			logFilename = arg.objectname + '.csv'
		else: 
			logFilename = spectrum.objectName + '.csv'
		if os.path.exists(logFilename):
			print "File exists!"
			logFile = open(logFilename, 'rt')
			for index, l in enumerate(logFile):
				if index != 0:
					oldLines.append(l)
			logFile.close()
			
		newLine = "%f, %f, %f\n"%(spectrum.HJD, velocity, velocityError)

		logLines.append(newLine)
		HJDString = "%f"%spectrum.HJD
		for l in oldLines:
			fields = l.split()[0]
			dateString =  fields.split(',')[0]
			
			if dateString == HJDString: 
				print "duplicate line, overwriting"
			else:
				logLines.append(l)
				
		logFile = open(logFilename, 'wt')
		logFile.write("HJD, velocity, velocity_error\n")

		logLines = sorted(logLines, reverse = False)
		print logLines
		for n in logLines:
			logFile.write(n)
			print n
			
		logFile.close()	
