#!/usr/bin/env python
import sys
import argparse
import subprocess
import ppgplot
import spectrumClasses
import scipy.optimize
import numpy
import timeClasses
import math

pathToCode = "/Users/rashley/astro/reductions/CSS081231/boris"
# pathToCode = "."

global iteration, mainPlotWindow, currentPlotWindow

def radians(angle):
	return angle/180. * math.pi

def replaceExpChar(value):
	if 'D' in value:
		value = value.replace('D', 'E')
	
	return value
	
def getModelSpectrum(angle, field, temperature, log_lambda, geometry):
	modelParams = {}
	modelParams['angle'] = float(angle)
	modelParams['field'] = field
	modelParams['temperature'] = temperature
	modelParams['log_lambda'] = log_lambda
	modelParams['geometry'] = geometry
	print "Creating a model for: ", modelParams
	
	filename = "test.dat"	
	parameterFile = open("ConstLambda_Ein", "w")
	parameterFile.write("Sichtwinkel      [grad] : %f\n"%angle)
	parameterFile.write("Magnetfeld       [MG]   : %f\n"%field)
	parameterFile.write("Temperatur       [keV]  : %f\n"%temperature)
	parameterFile.write("log(Lambda)      [1]    : %f\n"%log_lambda)
	parameterFile.write("Geometrie= 0 oder 1     : 0\n")
	parameterFile.close()

	modelCommand = [pathToCode + "/ConstLambda.bin"]
	outfile = open(filename, "w")
	subprocess.call(modelCommand, stdout = outfile)
	outfile.close()
					
	# Load the computed model
	print "Loading model:", filename		
	inputfile = open(filename, 'r')
	freqs = []
	i_0s = []
	i_1s = []
	i_s = []
	for line in inputfile:
		if line[0] == '#':
			continue
		if " warnung! " in line:
			continue
		data = line.split()
		replaceExpChar(data[0])
		freqs.append(float(replaceExpChar(data[0])))
		i_0s.append(float(replaceExpChar(data[1])))
		i_1s.append(float(replaceExpChar(data[2])))
		i_s.append(float(replaceExpChar(data[3])))		
			
		
	c = 3.E8
	angstroms = 1E-10
	wavelengths = [c / f for f in freqs]
	flambdas = []
	for fnu,l in zip(i_0s, wavelengths):
		flambda = fnu * c / (l*l)
		flambdas.append(flambda) 
		# print l, flambda
	wavelengths = [w / angstroms for w in wavelengths] 
	lowerWavelength = min(wavelengths)
	upperWavelength = max(wavelengths)
	lowerFlambda = min(flambdas)
	upperFlambda = max(flambdas)
	
	spectrum = spectrumClasses.spectrumObject()
	spectrum.setData(wavelengths, flambdas)
	spectrum.name = "model: angle=%f, temp=%f, field=%f"%(angle, temperature, field)
	spectrum.sortData()
	
	return spectrum
	
def getSampledModel(wavelengths, angle, field, temperature, log_lambda):
	global colour
	geometry = 0
	model = getModelSpectrum(angle, field, temperature, log_lambda, geometry)
	model.resample(wavelengths)
	modelArea = model.integrate()
	model.divide(modelArea)
	model.divide(1/observedArea)
	return model.flux
	
def quadratic(x, a0, a1, a2):
	y = a0*x*x + a1 *x + a2
	return y
	
def computeChiSq(spectrum, model):
	chi = 0
	for f1, f2 in zip(spectrum.flux, model):
		chi+= (f1-f2)**2 
	return chi
	
def getChiSqByParameters(params, *args):
	global iteration, mainPlotWindow, currentPlotWindow, colour, inclination, phase
	print "Params:", params
	beta = params[0]
	log_lambda = params[1]
	scale_factor = params[2]
	linear_offset = params[3]
	print "Args:", args
	temperature = args[0]
	field = args[1]
	
	# cos(theta) = cos(i)cos(beta) - sin(i)sin(beta)cos(phi + pi/2)
	cosTheta = math.cos(radians(inclination)) * math.cos(radians(beta)) - math.sin(radians(inclination)) * math.sin(radians(beta))*math.cos(phase + math.pi/2.)
	angle = math.acos(cosTheta) / math.pi * 180
	print "Angle: %f [deg], Field: %f [MG], Temperature:%f [keV], log_lambda: %f, scale: %f, offset: %f"%(angle, field, temperature, log_lambda, scale_factor, linear_offset)
	model = getSampledModel(observedSpectrum.wavelengths, angle, field, temperature, log_lambda)
	model = [m * scale_factor + linear_offset for m in model]
	chi = computeChiSq(observedSpectrum, model)
	allChiSqs.append(chi)
	print "Chi-squared:", chi
	startWavelength = min(observedSpectrum.wavelengths)
	endWavelength = max(observedSpectrum.wavelengths)
	
	# Draw the most recent iteration
	ppgplot.pgslct(currentPlotWindow)
	ppgplot.pgsci(1)
	ppgplot.pgenv(startWavelength, endWavelength, lowerFlux, upperFlux, 0, 0)
	ppgplot.pgline(observedSpectrum.wavelengths, observedSpectrum.flux)
	ppgplot.pgsci(4)
	ppgplot.pgline(observedSpectrum.wavelengths, model)
	ppgplot.pgsci(1)
	ppgplot.pglab("wavelength", "flux", "Current fit: %d"%iteration)
	
	# Overplot the iteration on the original diagram
	print "overplotting"
	ppgplot.pgslct(mainPlotWindow)
	ppgplot.pgsci(colour)
	ppgplot.pgline(observedSpectrum.wavelengths, model)
	colour += 1
	if colour>15: colour = 1
	ppgplot.pgsci(1)
	
	# Re-generate the Chi-Squared plot
	ppgplot.pgslct(chiSqPlotWindow)
	
	if iteration > 9:
		ppgplot.pgenv(0, iteration+1, 0, max(allChiSqs), 0, 0)
	else:
		ppgplot.pgenv(0, 10, 0, max(allChiSqs), 0, 0)
	iterations = range(iteration+1)
	ppgplot.pgpt(iterations, allChiSqs, 2)
	minCh = min(allChiSqs)
	medCh = numpy.median(allChiSqs)
	maxCh = max(allChiSqs)
	ppgplot.pglab("Iteration [n]", "Chi-squared", "Chi-squared values [%.2f, %.2f, %.2f]"%(minCh, medCh, maxCh))
	
	
	iteration += 1
	return chi
	
if __name__ == "__main__":
	iteration = 0
	parser = argparse.ArgumentParser(description='Fits a cyclotron model to an observed spectrum.')
	parser.add_argument('spectrum', type=str, help='JSON files containing the spectrum to be fit.')
	parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
	parser.add_argument('--angle', type=float, help='[Optional] Pole inclination (Beta) angle guess. Default is 60 degrees.')
	parser.add_argument('--temperature', type=float, help='[Optional] Plasma temperature in keV guess. Default is 24 keV.')
	parser.add_argument('--field', type=float, help='[Optional] Surface magnetic field strength MG guess. Default is 34 MG.')
	parser.add_argument('--loglambda', type=float,  help='[Optional] Log Lambda parameter guess. Default is 1.0.')
	parser.add_argument('-e', '--ephemeris', type=str, help = "Optional ephemeris file so we can calculate the phase of the spectrum.")
	 
	arg = parser.parse_args()
	print arg
	
	if arg.angle == None: angle = 80.   
	else: angle = arg.angle
	print "Magnetic pole inclination guess:", angle
	    
	if arg.field == None: field = 34.
	else: field = arg.field
	print "Magnetic field guess:", field
	
	if arg.temperature == None: temperature = 29.0
	else: temperature = arg.temperature
	print "Plasma temperature guess:", temperature
		
	if arg.loglambda == None: log_lambda = 1.0
	else: log_lambda = arg.loglambda
	print "Optical depth (log_lambda) guess:", log_lambda
	
	inclination = 81.3
	
	geometry = 0 			# Geometry 0 or 1
	colour = 1
		
	spectrum = spectrumClasses.spectrumObject()
	filename = arg.spectrum
	spectrum.loadFromJSON(filename)
	print "Loaded %s, contains %s."%(filename, spectrum.objectName)
	
	if arg.ephemeris!=None:
		HJD = spectrum.getProperty('HJD')
		print "HJD:", HJD
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.ephemeris)
		print ephemeris
		phase = ephemeris.getPhase(HJD)
		print "Phase:", phase
		phaseAngle = phase * 2 * math.pi 
	else: 
		print "We need an ephemeris to compute the phase angle. Exiting."
		sys.exit()
		
	phase = phaseAngle
		
	
	# Snip out Halpha 
	spectrum.snipWavelengthRange(6550, 6570)
	
	lowerWavelength = min(spectrum.wavelengths)
	upperWavelength = max(spectrum.wavelengths)
	observedSpectrumRange = (lowerWavelength, upperWavelength)
	lowerFlux = min(spectrum.flux)
	upperFlux = max(spectrum.flux)
	lowerFlux = 0
	
	mainPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgsci(1)
	ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
	ppgplot.pgline(spectrum.wavelengths, spectrum.flux)
	ppgplot.pglab("wavelength", "flux", spectrum.objectName)
	
	observedArea = spectrum.integrate()
	print "Wavelength range of observations:", observedSpectrumRange
		
	# modelPlotWindow = ppgplot.pgopen(arg.device)	
	# pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	# ppgplot.pgask(False)
	
	currentPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	ppgplot.pgslct(currentPlotWindow)
	ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
	ppgplot.pgline(spectrum.wavelengths, spectrum.flux)
	ppgplot.pglab("wavelength", "flux", "Current fit")
	
	chiSqPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	ppgplot.pgslct(chiSqPlotWindow)
	ppgplot.pgenv(0, 10, 0, 100, 0, 0)
	ppgplot.pglab("Iteration [n]", "Chi-squared", "Chi-squared values")
	
	
	colour = 2
	allChiSqs = []
	
	observedSpectrum = spectrum
	#         Angle      Temp     Scale factor  Linear offset
	bounds = [(60, 90), (10, 40), (0.5, 1.5),   (0.5, 1.5)]
	#        Field strength
	fixed = (36, 33)

	scale = 1.0
	offset = 0.1
	guess = [angle, log_lambda, scale, offset]
	fixed = (temperature, field)
	iteration = 0
	results = scipy.optimize.minimize(getChiSqByParameters, guess, args = fixed, method = 'Nelder-Mead')
	
	print results
	
