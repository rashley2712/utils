#!/usr/bin/env python
import sys
import argparse
import subprocess
import ppgplot
import spectrumClasses
import scipy.optimize

pathToCode = "/storage/astro2/phrnaw/reductions/CSS081231/boris"


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
		
	print "About to execute: " + str(modelCommand)
	subprocess.call(modelCommand)
				
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
	
def getSampledModel(wavelengths, angle, field, temperature, log_lambda, geometry):
	model = getModelSpectrum(angle, field, temperature, log_lambda, geometry)
	model.resample(wavelengths)
	modelArea = model.integrate()
	model.divide(modelArea)
	model.divide(1/observedArea)
	print "observed area", observedArea
	
	return model.flux

	
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
	
	modelPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	
	
	angle = 60.
	field = 34.
	temperature = 30.
	guess = [angle, field, temperature]
	
	modelSpectrum = getModelSpectrum(angle, field, temperature, 1, 0)
	lowerWavelength = min(modelSpectrum.wavelengths)
	upperWavelength = max(modelSpectrum.wavelengths)
	lowerFlambda = min(modelSpectrum.flux)
	upperFlambda = max(modelSpectrum.flux)
	lowerFlambda = 0

	ppgplot.pgslct(modelPlotWindow)
	ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlambda, upperFlambda, 0, 0)
	ppgplot.pglab("wavelength", "i_0", modelSpectrum.name)
		
	ppgplot.pgline(modelSpectrum.wavelengths, modelSpectrum.flux)

	
	#modelSpectrum.trimWavelengthRange(observedSpectrumRange[0], observedSpectrumRange[1])
	
	#modelArea = modelSpectrum.integrate()
	#modelSpectrum.divide(modelArea)
	#modelSpectrum.divide(1/observedArea)
	
	
	modelFlux = getSampledModel(spectrum.wavelengths, angle, field, temperature, 1, 0)
	
	ppgplot.pgslct(mainPGPlotWindow)
	# ppgplot.pgsci(colour)
	# ppgplot.pgline(modelSpectrum.wavelengths, modelSpectrum.flux)
	# colour+= 1
	ppgplot.pgsci(colour)
	ppgplot.pgline(spectrum.wavelengths, modelFlux)
	colour+= 1
	if colour > 15: colour = 1 
	

	
