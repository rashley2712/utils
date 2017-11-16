#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import loadingSavingUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot
import gSheets


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
	# w = 3.4
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

class dataLog:
	def __init__(self):
		self.measurements = []

	def addMeasurement(self, HJD, velocity, velocityError, fwhm, wavelength):
		measurement = {}
		measurement['HJD'] = HJD
		measurement['velocity'] = velocity
		measurement['velocityError'] = velocityError
		measurement['fwhm'] = fwhm
		measurement['wavelength'] = wavelength
		
		# Check if the measurement is a duplicate (ie same HJD)
		found = -1
		for index, m in enumerate(self.measurements):
			if m['HJD'] == HJD: 
				found = index
				m['velocity'] = velocity
				m['velocityError'] = velocityError
				m['fwhm'] = fwhm
				m['wavelength'] = wavelength
		if found == -1:
			self.measurements.append(measurement)
	
		
	def writeToFile(self, filename):
		logFile = open(filename, 'wt')
		logFile.write("HJD, velocity, velocity_error, fwhm, wavelength\n")
		for m in self.measurements: 
			outString = "%10.10f, %10.10f, %10.10f, %10.10f, %10.10f\n"%(m['HJD'], m['velocity'], m['velocityError'], m['fwhm'], m['wavelength'])
			logFile.write(outString)
		logFile.close()	
		
	def getSavedValues(self, HJD):
		for m in self.measurements:
			if (float(m['HJD']) == float(HJD)): 
				return (m['wavelength'], m['fwhm'])
		
		return (-1, -1)
		
	def sortByHJD(self):
		self.measurements =  sorted(self.measurements, key=lambda object: object['HJD'], reverse = False)
		
	def loadFromFile(self, filename):
		if not os.path.exists(filename): return
		inputFile = open(filename, 'rt')
		for line in inputFile:
			parts = line.strip().split(',')
			if parts[0].strip(',') == 'HJD': continue 
			HJD = float(parts[0].strip(','))
			velocity = float(parts[1].strip(','))
			velocityError = float(parts[2].strip(','))
			fwhm = float(parts[3].strip(','))
			wavelength = float(parts[4].strip(','))
			self.addMeasurement(HJD, velocity, velocityError, fwhm, wavelength)
		inputFile.close()
		
	def __str__(self):
		if len(self.measurements)==0: return
		retStr = "HJD, velocity, velocityErr, fwhm, wavelength\n"
		for m in self.measurements: 
			retStr+="%10.10f, %f, %f, %f, %f\n"%(m['HJD'], m['velocity'], m['velocityError'], m['fwhm'], m['wavelength'])
		return retStr

if __name__ == "__main__":
	NA_labwavelength = 8183.2556
	
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
	parser.add_argument('--fixwidth', action='store_true', help='Use the width value that is stored in the Google sheet. ')
	parser.add_argument('--skipgood', action='store_true', help='Don''t try to fit spectra that are marked as ''good'' in the sheet.')
	parser.add_argument('-o', '--objectname', type=str, help='[Optional] Object name for the output log file.')
	arg = parser.parse_args()
	
	# docsCredentials = gSheets.get_credentials()
	# sampleData = gSheets.getSampleData("1BxiMVs0XRA5nFMdKvBdBZjgmUUqptlbs74OgvE2upms")
	
	docInstance = gSheets.gSheetObject()
	docInstance.initCredentials()
	docInstance.setDocID('11fsbzSII1u1-O6qQUB8P0RzvJ8MzC5VHIASsZTYplXc')
	docInstance.setObjectName(arg.objectname)
	docInstance.loadAllReadings()
	# docInstance.createSheet("testname")
	
	
	# print arg
	defaultWidth = 1.0 # Default width of line in Angstrom
	defaultWavelength = NA_labwavelength # Blueward doublet line in Angstrom
	defaultVelocity = 0.0
	defaultVelocityError = 0.0 
	
	
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
	
	
	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
		
	"""# Load any existing data from the dataLog file
	recordedData = dataLog()
	if arg.objectname is not None:
		logFilename = arg.objectname + '.csv'
		recordedData.loadFromFile(logFilename)
		
	recordedData.sortByHJD()
	# print recordedData
	"""
	mainPGPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	yUpper = 2.5
	yLower = -0.5
	
	fitPGPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	
	for spectrum in spectra:
		ppgplot.pgslct(mainPGPlotWindow)
		
		ppgplot.pgsci(1)
		lowerWavelength = min(spectrum.wavelengths)
		upperWavelength = max(spectrum.wavelengths)
		lowerFlux = min(spectrum.flux)
		upperFlux = max(spectrum.flux)
		lowerLimit = min(spectrum.fluxErrors)
		ppgplot.pgenv(lowerWavelength, upperWavelength, 0, upperFlux, 0, 0)
		ppgplot.pgbin(spectrum.wavelengths, spectrum.flux)
		ppgplot.pgbin(spectrum.wavelengths, spectrum.fluxErrors)
		ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "flux [%s]"%spectrum.fluxUnits, "%s [%s]"%(spectrum.objectName, spectrum.loadedFromFilename))
		
		# Grab the continuum from either side of the spectrum
		
		lowerCut = 8175
		upperCut = 8205
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
		y_errors = continuumSpectrum.fluxErrors
		results, covariance = scipy.optimize.curve_fit(quad, x_values, y_values, guess, )
		errors = numpy.sqrt(numpy.diag(covariance))
		# print "quadratic result:", results
		# print "quadratic errors:", errors
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
			normalisedSpectrum.fluxErrors[index] = spectrum.fluxErrors[index]/yFit[index]
		
		lowerWavelength = min(normalisedSpectrum.wavelengths)
		upperWavelength = max(normalisedSpectrum.wavelengths)
		lowerFlux = min(normalisedSpectrum.flux)
		upperFlux = max(normalisedSpectrum.flux)
	
		ppgplot.pgslct(fitPGPlotWindow)
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
		ppgplot.pgsci(1)
		ppgplot.pgsls(3)
		ppgplot.pgbin(normalisedSpectrum.wavelengths, normalisedSpectrum.flux)
		ppgplot.pgsls(1)
		ppgplot.pgsci(2)
		ppgplot.pgbin(normalisedSpectrum.wavelengths, normalisedSpectrum.fluxErrors + lowerFlux)
		ppgplot.pgsci(1)
		ppgplot.pglab("wavelength [%s]"%spectrum.wavelengthUnits, "flux [normalised]", "%s [%s]"%(spectrum.objectName, spectrum.loadedFromFilename))
		
		featureSpectrum = copy.deepcopy(normalisedSpectrum)
		featureSpectrum.trimWavelengthRange(lowerCut, upperCut)
		ppgplot.pgbin(featureSpectrum.wavelengths, featureSpectrum.flux)
		
		"""# Check if there is a guess value to use for the wavelength and the width of the fit
		wavelength, fwhm = recordedData.getSavedValues(spectrum.HJD)
		if (wavelength!=-1):
			centroidWavelength = wavelength
			width = fwhm
			print "Found previously saved values for %10.10f: wavelength = %fAA and fwhm = %fAA"%(spectrum.HJD, centroidWavelength, fwhm)
			constant = 1.0
			depth = -0.5
			"""
		if (docInstance.hasReadingFor(spectrum.HJD)):
			print "Found a previous fit for this spectrum %f"%spectrum.HJD 
			(width, wavelength) = docInstance.getFitByHJD(spectrum.HJD)
			if arg.skipgood and docInstance.getGoodFlag(spectrum.HJD):
				print "Skipping spectrum as it is marked as a good fit"
				continue
			a0 = 1.0
			a1 = -0.5
			a2 = wavelength
			a3 = width
			print "Using width: %f and wavelength: %f"%(width, wavelength)
			constant = a0
			depth = a1
			centroidWavelength = a2
			fixWidth = width
			
		else: 
			print "No previous fit found for this spectrum %f"%spectrum.HJD 
			docInstance.addNewMeasurement(spectrum.HJD, defaultVelocity, defaultVelocityError, defaultWidth, defaultWavelength, False)
			docInstance.writeAllReadings()
			# Fit a single gaussian to the Na doublet blue line 
			a0 = 1.0    	# Constant term  
			a1 = -0.5		# 'Depth' of the line
			a2 = defaultWavelength		# Wavelength of the centre of the line
			a3 = defaultWidth		# Width of the line
		
		
			guess = numpy.array([a0, a1, a2, a3])
			x_values = featureSpectrum.wavelengths
			y_values = featureSpectrum.flux
			y_errors = featureSpectrum.fluxErrors
			results, covariance = scipy.optimize.curve_fit(gaussian, x_values, y_values, guess, y_errors, absolute_sigma = True)
			errors = numpy.sqrt(numpy.diag(covariance))
			a0 = results[0]
			a1 = results[1]
			a2 = results[2]
			a3 = results[3]
			print "Centroid wavelength %f [%f]"%(a2, errors[2])
			print "Width %f [%f]"%(a3, errors[3])
			xFit = spectrum.wavelengths
			yFit = gaussian(numpy.array(xFit), a0, a1, a2, a3)
		
		width = a3
		if arg.fixwidth:
			print "Using the width as specified in the sheet ... %f angstrom"%fixWidth
			width = fixWidth
		centroidWavelength = a2
		depth = a1
		constant = a0

		# Draw the first single Gaussian fit
		#currentColour = ppgplot.pgqci()
		#ppgplot.pgsci(3)
		#ppgplot.pgline(xFit, yFit)
		#ppgplot.pgsci(currentColour)

		# Now fit the double gaussian with a fixed width and separation
		a0 = constant    			# Constant term  
		a1 = depth					# 'Depth' of the line
		a2 = centroidWavelength		# Wavelength of the centre of the blueward line
		guess = numpy.array([a0, a1, a2])
		x_values = featureSpectrum.wavelengths
		y_values = featureSpectrum.flux
		y_errors = [e for e in featureSpectrum.fluxErrors]
		#print "Fluxes:", y_values
		#print "Flux errors:", y_errors
		results, covariance = scipy.optimize.curve_fit(doubleGaussian, x_values, y_values, guess, y_errors, absolute_sigma = True)
		errors = numpy.sqrt(numpy.diag(covariance))
		# print "double gaussian result:", results
		# print "double gaussian errors:", errors
		a0 = results[0]
		a1 = results[1]
		a2 = results[2]
		wavelength = a2
		wavelengthError = errors[2]
		print "Centroid blueward wavelength %f [%f]"%(wavelength, wavelengthError)
		velocity = (wavelength - NA_labwavelength)/NA_labwavelength * 3E5 
		velocityError = 3E5 /NA_labwavelength * wavelengthError
		print "Velocity %f [%f]"%(velocity, velocityError)

		xFit = spectrum.wavelengths
		yFit = doubleGaussian(numpy.array(xFit), a0, a1, a2)
		ppgplot.pgsci(3)
		ppgplot.pgline(xFit, yFit)

		# Compute ChiSquared of the fit
		chiSq = 0
		sigma = 0
		for x, flux, fluxError in zip(x_values, y_values, y_errors):
			fittedFlux = doubleGaussian(x, a0, a1, a2) 			
			# print wavelength, flux, fittedFlux, fluxError
			chiSq+= ((flux - fittedFlux)/fluxError)**2
			sigma+= (flux - fittedFlux)**2
		print "Chi squared", chiSq
		sigma = numpy.sqrt(  sigma/(len(x_values) - 1)    )
		reducedChiSq = chiSq / (len(x_values) - 3)
		print "sigma:", sigma
		print "Reduced Chi squared", reducedChiSq

		threeSigmaPlus = [f + 3*sigma for f in yFit]
		threeSigmaMinus = [f - 3*sigma for f in yFit]

		newX = []
		newY = []
		newYE = []
		redo = False
		for w, f, fe in zip(x_values, y_values, y_errors):
			fittedFlux = doubleGaussian(w, a0, a1, a2)
			if abs(f - fittedFlux) > 3*sigma: 
				print w, f, 'is more than 3 sigma from the fit', fittedFlux, 'rejecting it'
				continue
			else:
				newX.append(w)
				newY.append(f)
				newYE.append(fe)
				redo = True 

		print len(x_values)
		x_values = newX
		print len(x_values)
		y_values = newY
		y_errors = newYE
		if redo:
			print "Redoing the fit with the leftover points"
			guess = numpy.array([a0, a1, a2])
			results, covariance = scipy.optimize.curve_fit(doubleGaussian, x_values, y_values, guess, y_errors, absolute_sigma = True)	
			errors = numpy.sqrt(numpy.diag(covariance))
			a0 = results[0]
			a1 = results[1]
			a2 = results[2]
			wavelength = a2
			wavelengthError = errors[2]
			print "Centroid blueward wavelength %f [%f]"%(wavelength, wavelengthError)
			velocity = (wavelength - NA_labwavelength)/NA_labwavelength * 3E5 
			velocityError = 3E5 /NA_labwavelength * wavelengthError
			print "Velocity %f [%f]"%(velocity, velocityError)
	
		
		ppgplot.pgsls(2)
		ppgplot.pgline(xFit, threeSigmaPlus)
		ppgplot.pgline(xFit, threeSigmaMinus)

		ppgplot.pgsci(5)
		ppgplot.pgsls(2)
		ppgplot.pgline([NA_labwavelength, NA_labwavelength], [lowerFlux, upperFlux])
		ppgplot.pgsci(6)
		ppgplot.pgline([lowerCut, lowerCut], [lowerFlux, upperFlux])
		ppgplot.pgline([upperCut, upperCut], [lowerFlux, upperFlux])
		ppgplot.pgsls(1)
		ppgplot.pgsci(1)
		
		if numpy.isinf(velocityError): velocityError = 9E9
		if not query_yes_no("Are you happy with the fit?"):
			print "Saving the value with the flag raised."
			docInstance.addNewMeasurement(spectrum.HJD, velocity, velocityError, width, wavelength, False)
			docInstance.writeAllReadings()
		else:
			docInstance.addNewMeasurement(spectrum.HJD, velocity, velocityError, width, wavelength, True)
			docInstance.writeAllReadings()
		"""recordedData.addMeasurement(spectrum.HJD, velocity, velocityError, width, wavelength)
		recordedData.sortByHJD()
		recordedData.writeToFile(logFilename)
		"""
		
	# Write all to a CSV file
	data = docInstance.readings
	outputLog = open(arg.objectname + ".csv", "wt")
	outputLog.write("HJD, Velocity, VelErr, Width, Wavelength, Good\n")
	for d in data:
		print d
		if d['good'] == 1:
			outputLog.write("%f, %f, %f, %f, %f, %f\n"%(d['HJD'], d['RV'], d['RV error'], d['width'], d['wavelength'], d['good']))
	outputLog.close()
		
