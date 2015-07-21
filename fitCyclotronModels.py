#!/usr/bin/env python
import sys
import argparse
import subprocess
import ppgplot
import spectrumClasses

def replaceExpChar(value):
	if 'D' in value:
		value = value.replace('D', 'E')
	
	return value
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Interactively fits a cyclotron model to an observed spectrum.')
	parser.add_argument('spectrum', type=str, help='JSON files containing the spectrum to be fit.')
	parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
	 
	arg = parser.parse_args()
	print arg
	
	pathToCode = "/storage/astro2/phrnaw/reductions/CSS081231/boris"
	
	angle = 60.0    		# Sight angle in degrees
	field = 64.0   			# Field strength in MG
	temperature =20.0		# Temperature in keV
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
	
	for angle in range(90, 10, -10):
		modelParams = {}
		modelParams['angle'] = float(angle)
		modelParams['field'] = field
		modelParams['temperature'] = temperature
		modelParams['log_lambda'] = log_lambda
		modelParams['geometry'] = geometry
		print modelParams
		
		filename = "test"	
		modelCommand = [pathToCode + "/ConstLambda"]
		modelCommand.append(str(field))
		modelCommand.append(str(angle))
		modelCommand.append(str(temperature))
		modelCommand.append(str(log_lambda))
		modelCommand.append(str(geometry))
		modelCommand.append(str(filename))
		
		print "About to execute: " + str(modelCommand)
		subprocess.call(modelCommand)
				
		# Load and plot the model
		filename+= ".dat"
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
			
			
		lowerFreq = min(freqs)
		upperFreq = max(freqs)
		upper_i0 = max(i_0s)
		lower_i0 = min(i_0s)
		upper_i1 = max(i_1s)
		lower_i1 = min(i_1s)
		upper_i = max(i_s)
		lower_i = min(i_s)
		
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
			
		# trimLower = 5000
		# trimUpper = 10000
		# newWavelengths = []
		# newFlambdas = []
		# for f, l in zip(flambdas, wavelengths):
		#	if l > trimLower and l < trimUpper:
		#		newWavelengths.append(l)
		#		newFlambdas.append(f)
					
		# Now create a spectrum object from our model
		modelSpectrum = None
		modelSpectrum = spectrumClasses.spectrumObject()
		modelSpectrum.setData(wavelengths, flambdas)
		modelSpectrum.name = "model: angle %f"%angle
		modelSpectrum.sortData()
		
		modelSpectrum.trimWavelengthRange(observedSpectrumRange[0], observedSpectrumRange[1])
		
		modelArea = modelSpectrum.integrate()
		modelSpectrum.divide(modelArea)
		modelSpectrum.divide(1/observedArea)
		
		
		
		lowerWavelength = min(modelSpectrum.wavelengths)
		upperWavelength = max(modelSpectrum.wavelengths)
		lowerFlambda = min(modelSpectrum.flux)
		upperFlambda = max(modelSpectrum.flux)
		lowerFlambda = 0
		
		ppgplot.pgslct(modelPlotWindow)
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlambda, upperFlambda, 0, 0)
		ppgplot.pglab("wavelength", "i_0", modelSpectrum.name)
			
		ppgplot.pgline(modelSpectrum.wavelengths, modelSpectrum.flux)
		
		ppgplot.pgslct(mainPGPlotWindow)
		ppgplot.pgsci(2)
		ppgplot.pgline(modelSpectrum.wavelengths, modelSpectrum.flux)
		
	
	
