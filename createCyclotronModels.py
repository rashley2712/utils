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

	parser = argparse.ArgumentParser(description='Creates a series of cyclotron models using the ConstLambda code.')
	parser.add_argument('--device', type=str, help='[Additional plot device] PGPLOT device. "/xs" is always selected.')
	parser.add_argument('-t', '--temperature', type=float, default='10.0', help = 'Temperature in keV. Default is 10 keV.')
	parser.add_argument('-b', '--field', type=float, default='34.0', help = 'Surface field strength in megagauss (MG). Default is 34 MG.')
	 
	arg = parser.parse_args()
	print arg
	
	xDevice = "/xs"
	
	pathToCode = "/storage/astro2/phrnaw/reductions/CSS081231/boris"
	
	angle = 90    					# Sight angle in degrees
	field = arg.field     			# Field strength in MG
	temperature = arg.temperature	# Temperature in keV
	log_lambda = 1      			# log(lambda)
	geometry = 0 					# Geometry 0 or 1
	
	device = arg.device
	# device = "models.ps/ps"
	
	mainPGPlotWindow = ppgplot.pgopen(xDevice)	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgask(False)
		
	modelSpectra = []
	
	for angle in range(10, 95, 5):
		# filename = "B%fW%fT%fL%f.dat"%(angle, field, temperature, log_lambda)
		filename = "test"
		
		parameterFile = open("ConstLambda_Ein", "w")
		parameterFile.write("Sichtwinkel      [grad] : %f\n"%angle)
		parameterFile.write("Magnetfeld       [MG]   : %f\n"%field)
		parameterFile.write("Temperatur       [keV]  : %f\n"%temperature)
		parameterFile.write("log(Lambda)      [1]    : %f\n"%log_lambda)
		parameterFile.write("Geometrie= 0 oder 1     : 0\n")
		parameterFile.close()
		modelCommand = [pathToCode + "/ConstLambda.bin"]
		
		print "About to execute: " + str(modelCommand)
		outfile = open("test.dat", "w")
		subprocess.call(modelCommand, stdout = outfile)
		outfile.close()
		
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
		
		modelSpectrum = spectrumClasses.spectrumObject()
		modelSpectrum.setData(wavelengths, flambdas)
		modelSpectrum.trimWavelengthRange(5000, 9000)
		modelSpectrum.angle = angle
		modelSpectra.append(modelSpectrum)
		
		lowerWavelength = min(modelSpectrum.wavelengths)
		upperWavelength = max(modelSpectrum.wavelengths)
		lowerFlambda = min(modelSpectrum.flux)
		upperFlambda = max(modelSpectrum.flux)
		lowerFlambda = 0
		
		lowerFlux = 0
		upperFlux = 0
		
		for s in modelSpectra:
			fmax = max(s.flux)
			if fmax>upperFlux:
				upperFlux = fmax
		
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
		colour = 1
		for s in modelSpectra:
			ppgplot.pgsci(colour)
			ppgplot.pgline(s.wavelengths, s.flux)
			plotx = 8500
			ploty = s.getNearestFlux(plotx)
			ppgplot.pgptxt(plotx, ploty, 0, 0, "%2.0f'"%(s.angle)) 
			colour+= 1
			if colour>15: colour = 1
			
		ppgplot.pgsci(1)	
		label = "B=%.0f MG, T=%.1f keV"%(field, temperature)
		ppgplot.pglab("wavelength", "i_0", label)
		
	ppgplot.pgclos()

	if arg.device!=None:
		# Write to an additional output devide	
		mainPGPlotWindow = ppgplot.pgopen(arg.device)	
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pgask(False)
		ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
		colour = 1
		for s in modelSpectra:
			ppgplot.pgsci(colour)
			ppgplot.pgline(s.wavelengths, s.flux)
			plotx = 8500
			ploty = s.getNearestFlux(plotx)
			ppgplot.pgptxt(plotx, ploty, 0, 0, "%2.0f'"%(s.angle)) 
			colour+= 1
			if colour>15: colour = 1
			
		ppgplot.pgsci(1)	
		label = "B=%.0f MG, T=%.1f keV"%(field, temperature)
		ppgplot.pglab("wavelength", "i_0", label)
		ppgplot.pgclos()
