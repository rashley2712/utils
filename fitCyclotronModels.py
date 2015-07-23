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
	
def getSampledModel(wavelengths, angle, field, temperature):
	log_lambda = 1.0
	geometry = 0
	model = getModelSpectrum(angle, field, temperature, log_lambda, geometry)
	model.resample(wavelengths)
	modelArea = model.integrate()
	model.divide(modelArea)
	model.divide(1/observedArea)
	ppgplot.pgsci(colour)
	ppgplot.pgline(model.wavelengths, model.flux)
	print "colour:", colour
	colour+- 1
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
	angle = params[0]
	temperature = params[1]
	scale_factor = params[2]
	linear_offset = params[3]
	field = args[0]
	print "Angle: %f [deg], Field: %f [MG], Temperature:%f [keV], scale: %f, offset: %f"%(angle, field, temperature, scale_factor, linear_offset)
	model = getSampledModel(observedSpectrum.wavelengths, angle, field, temperature)
	model = [m * scale_factor + linear_offset for m in model]
	chi = computeChiSq(observedSpectrum, model)
	print "Chi-squared:", chi
	return chi
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Interactively fits a cyclotron model to an observed spectrum.')
	parser.add_argument('spectrum', type=str, help='JSON files containing the spectrum to be fit.')
	parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
	 
	arg = parser.parse_args()
	print arg
	
	
	angle = 60.0    		# Sight angle in degrees
	field = 34.0   			# Field strength in MG
	temperature =35.0		# Temperature in keV
	log_lambda = 1      	# log(lambda)
	geometry = 0 			# Geometry 0 or 1
	
		
	spectrum = spectrumClasses.spectrumObject()
	filename = arg.spectrum
	spectrum.loadFromJSON(filename)
	print "Loaded %s, contains %s."%(filename, spectrum.objectName)
		
	mainPGPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(True)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	yUpper = 2.5
	yLower = -0.5
	ppgplot.pgsci(1)
	lowerWavelength = min(spectrum.wavelengths)
	upperWavelength = max(spectrum.wavelengths)
	observedSpectrumRange = (lowerWavelength, upperWavelength)
	lowerFlux = min(spectrum.flux)
	upperFlux = max(spectrum.flux)
	lowerFlux = 0
	ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
	ppgplot.pgline(spectrum.wavelengths, spectrum.flux)
	ppgplot.pglab("wavelength", "flux", spectrum.objectName)
	
	observedArea = spectrum.integrate()
	print "Wavelength range of observations:", observedSpectrumRange
	
	"""spectrum.trimWavelengthRange(6000, 6900)
	a0 = -2.615E-6
	a1 = 3.373E-2
	a2 = -1.08E2
	yFit = quadratic(numpy.array(spectrum.wavelengths), a0, a1, a2)
	ppgplot.pgline(spectrum.wavelengths, yFit)
	
	x = numpy.array(spectrum.wavelengths)
	da0 = .1E-6
	da1 = .1E-2
	da2 = .1E2
	
	chisq
	alpha00 = (quadratic(x, a0 + da0, a1, a2) - 2 * quadratic(x, a0, a1, a2) + quadratic(x, a0 - da0, a1, a2)) / (da0*da0)
	alpha11 = (quadratic(x, a0, a1 + da1, a2) - 2 * quadratic(x, a0, a1, a2) + quadratic(x, a0, a1 - da1, a2)) / (da1*da1)
	alpha22 = (quadratic(x, a0, a1, a2 + da2) - 2 * quadratic(x, a0, a1, a2) + quadratic(x, a0, a1, a2 - da2)) / (da2*da2)
	alpha01 = (quadratic(x, a0 + da0, a1 + da1, a2) - quadratic(x, a0 - da0, a1 + da1, a2) - quadratic(x, a0+da0, a1 - da1, a2) + quadratic(x, a0-da0, a1-da1, a2)) / (4*da0*da1)
	alpha = [ [alpha00, alpha01, 0], [alpha01, alpha22, 0], [0, 0, alpha22] ]
	print alpha
	
	
	sys.exit()
	"""
	
	modelPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	
	
	angle = 60.
	field = 34.
	temperature = 20.
	guess = [angle, field, temperature]
	
	modelSpectrum = getModelSpectrum(angle, field, temperature, 1, 0)
	lowerWavelength = min(modelSpectrum.wavelengths)
	upperWavelength = max(modelSpectrum.wavelengths)
	lowerFlambda = min(modelSpectrum.flux)
	upperFlambda = max(modelSpectrum.flux)
	lowerFlambda = 0

	colour = 1
	ppgplot.pgslct(modelPlotWindow)
	ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlambda, upperFlambda, 0, 0)
	ppgplot.pglab("wavelength", "i_0", modelSpectrum.name)
		
	ppgplot.pgline(modelSpectrum.wavelengths, modelSpectrum.flux)
	
	
	ppgplot.pgslct(mainPGPlotWindow)
	ppgplot.pgsci(2)
	
	# y_errors = numpy.ones(len(spectrum.flux))
	colour = 3
	
	"""for i in range(10):
		print "Iteration:", i
		model = getSampledModel(spectrum.wavelengths, angle, field, temperature)
		chi_sq = computeChiSq(spectrum, model)
		dt = 1
		chi_sqplusdt = computeChiSq(spectrum, getSampledModel(spectrum.wavelengths, angle, field, temperature + dt))
		chi_sqminusdt = computeChiSq(spectrum, getSampledModel(spectrum.wavelengths, angle, field, temperature - dt))
		dXdt = (chi_sqplusdt - chi_sqminusdt ) / (2*dt)

		da = 1
		chi_sqplusda = computeChiSq(spectrum, getSampledModel(spectrum.wavelengths, angle + da, field, temperature))
		chi_sqminusda = computeChiSq(spectrum, getSampledModel(spectrum.wavelengths, angle + da, field, temperature))
		dXda = (chi_sqplusda - chi_sqminusda ) / (2*da)


		print "ChiSq", chi_sq
		print "dChiSq/dt", dXdt
	"""	
	
	observedSpectrum = spectrum
	#         Angle      Temp     Scale factor  Linear offset
	bounds = [(45, 80), (10, 40), (0.5, 1.5),   (0.5, 1.5)]
	#        Field strength
	fixed = (36, 33)

	results = scipy.optimize.differential_evolution(getChiSqByParameters, bounds, args = fixed)
	print results
	
