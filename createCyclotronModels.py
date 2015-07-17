#!/usr/bin/env python
import sys
import argparse
import subprocess
import ppgplot

def replaceExpChar(value):
	if 'D' in value:
		value = value.replace('D', 'E')
	
	return value
	
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Creates a series of cyclotron models using the ConstLambda code.')
	 
	arg = parser.parse_args()
	print arg
	
	pathToCode = "/storage/astro2/phrnaw/reductions/CSS081231/boris"
	
	angle = 90    		# Sight angle in degrees
	field = 34    		# Field strength in MG
	temperature =50		# Temperature in keV
	log_lambda = 1      # log(lambda)
	geometry = 0 		# Geometry 0 or 1
	
	device = "/xs"
	# device = "models.ps/ps"
	
	mainPGPlotWindow = ppgplot.pgopen(device)	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgask(False)
	setEnv = False
	colour = 1
		
	for angle in range(90, 0, -5):
		print field, angle
		# filename = "B%fW%fT%fL%f.dat"%(angle, field, temperature, log_lambda)
		filename = "test"
		
		modelCommand = [pathToCode + "/ConstLambda"]
		modelCommand.append(str(angle))
		modelCommand.append(str(field))
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
		lower_i0 = max(i_0s)
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
		wavelengths = [w / angstroms for w in wavelengths] 
		lowerWavelength = min(wavelengths)
		upperWavelength = max(wavelengths)
		lowerFlambda = min(flambdas)
		upperFlambda = max(flambdas)
		
		trimLower = 5000
		trimUpper = 10000
		newWavelengths = []
		newFlambdas = []
		for f, l in zip(flambdas, wavelengths):
			if l > trimLower and l < trimUpper:
				newWavelengths.append(l)
				newFlambdas.append(f)
		
		wavelengths = newWavelengths
		flambdas = newFlambdas
		
		lowerWavelength = min(wavelengths)
		upperWavelength = max(wavelengths)
		lowerFlambda = min(flambdas)
		upperFlambda = max(flambdas)
		lowerFlambda = 0
				
		if not setEnv:
			ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlambda, upperFlambda, 0, 0)
			ppgplot.pglab("wavelength", "i_0", "")
			setEnv = True
		ppgplot.pgsci(colour)
		colour+= 1
		
		ppgplot.pgline(wavelengths, flambdas)
			
			
