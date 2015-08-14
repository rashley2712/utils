#!/usr/bin/env python
import sys
import argparse
import subprocess
import ppgplot
import spectrumClasses
import scipy.optimize
import numpy

pathToCode = "/storage/astro2/phrnaw/reductions/CSS081231/boris"
pathToCode = "."

global iteration, mainPlotWindow, currentPlotWindow

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
	global iteration, mainPlotWindow, currentPlotWindow, colour
	angle = params[0]
	temperature = params[1]
	scale_factor = params[2]
	linear_offset = params[3]
	field = params[4]
	log_lambda = params[5]
	print "Angle: %f [deg], Field: %f [MG], Temperature:%f [keV], log_lambda: %f, scale: %f, offset: %f"%(angle, field, temperature, log_lambda, scale_factor, linear_offset)
	model = getSampledModel(observedSpectrum.wavelengths, angle, field, temperature, log_lambda)
	model = [m * scale_factor + linear_offset for m in model]
	chi = computeChiSq(observedSpectrum, model)
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
	
	
	iteration += 1
	return chi
	
if __name__ == "__main__":
	iteration = 0
	parser = argparse.ArgumentParser(description='Fits a cyclotron model to an observed spectrum.')
	parser.add_argument('spectrum', type=str, help='JSON files containing the spectrum to be fit.')
	parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
	parser.add_argument('--angle', type=float, default = 60.0, help='[Optional] Line of sight angle guess. Default is 60 degrees.')
	parser.add_argument('--temperature', type=float, default = 24.0, help='[Optional] Plasma temperature in keV guess. Default is 24 keV.')
	parser.add_argument('--field', type=float, default = 34.0, help='[Optional] Surface magnetic field strength MG guess. Default is 34 MG.')
	parser.add_argument('--loglambda', type=float, default = 1.0, help='[Optional] Log Lambda parameter guess. Default is 1.0.')
	 
	arg = parser.parse_args()
	print arg
	
	angle = 10.0    		# Sight angle in degrees
	field = 34.0   			# Field strength in MG
	temperature =35.0		# Temperature in keV
	log_lambda = 1      	# log(lambda)
	geometry = 0 			# Geometry 0 or 1
	colour = 1
		
	spectrum = spectrumClasses.spectrumObject()
	filename = arg.spectrum
	spectrum.loadFromJSON(filename)
	print "Loaded %s, contains %s."%(filename, spectrum.objectName)
	
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
	
	angle = 60.
	field = 34.
	temperature = 20.
	guess = [angle, field, temperature]
	
	
	colour = 2
	
	observedSpectrum = spectrum
	#         Angle      Temp     Scale factor  Linear offset
	bounds = [(45, 80), (10, 40), (0.5, 1.5),   (0.5, 1.5)]
	#        Field strength
	fixed = (36, 33)

	angle = arg.angle
	field = arg.field
	temperature = arg.temperature
	log_lambda = arg.loglambda
	guess = [angle, temperature, 1.0, 0.1, field, log_lambda]
	iteration = 0
	results = scipy.optimize.minimize(getChiSqByParameters, guess, args = fixed, method = 'Nelder-Mead')
	
	print results
	
