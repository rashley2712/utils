#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import numpy, math
import matplotlib.pyplot
import argparse
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy

def quad(x, a1, a2, a3):
	y = a1 * x * x + a2 *x + a3
	return y

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a set of model spectra from the Sloane model set and plots them.')
	parser.add_argument('inputFiles', type=str, nargs='+', help='Molly files containing the spectra.')
	parser.add_argument('--list', action = 'store_true', help='The input file is a list of filenames.')
	parser.add_argument('-s', '--spectrum', type=str, help='Optional spectrum file to compare to the models.')
		 
	arg = parser.parse_args()
	print arg
	
	if arg.list:
		fileList = []
		print "Opening", arg.inputFiles[0]
		readfileList = open(arg.inputFiles[0], 'r')
		for line in readfileList:
			fileList.append(line.strip('\n'))  
	else:
		fileList = arg.inputFiles
	spectra = []

	for f in fileList:
		print "Loading", f
		spectrum = spectrumClasses.spectrumObject()
		spectrum.objectName = f
		inputFile = open(f, 'r')
		wavelengths = []
		fluxes = []
		
		for line in inputFile:
			if line[0] == '#':
				continue;
			params =  line.split()
			wavelength = float(params[0])
			flux = float(params[1])
			wavelengths.append(wavelength)
			fluxes.append(flux)

		spectrum.wavelengths = wavelengths
		spectrum.flux = fluxes
		spectra.append(spectrum)
		inputFile.close()

	
	numSpectra = len(spectra)
	
	if arg.spectrum!=None:
		columnNames, data = loadingSavingUtils.loadNewCSV(arg.spectrum)

		observedSpectrum = spectrumClasses.spectrumObject()
		observedSpectrum.objectName = "CSS081231"
		observedSpectrum.setData(data['wavelength'], data['flux'])

		tioSpectrum = copy.deepcopy(observedSpectrum)
		tioSpectrum.trimWavelengthRange(7000, 7500)
		tioArea = tioSpectrum.integrate()
		# observedSpectrum.divide(tioArea)
		lower, upper = observedSpectrum.wavelengthRange
	


	lowest_chi = 100;
	bestFit = "none"
	bestIndex = 0
	matplotlib.pyplot.figure(figsize=(20, 4*1 + 1))
	for index, s in enumerate(spectra):
		matplotlib.pyplot.clf()
		s.resample(observedSpectrum.wavelengths)
		subSpectrum = copy.deepcopy(s)
		newLength = subSpectrum.trimWavelengthRange(7000, 7500)
		area = subSpectrum.integrate()
		s.divide(area)
		s.divide(1/tioArea)
		# s.trimWavelengthRange(lower, upper)
		# Compute chi=squared
		chi_sq = 0
		for (data_w, data_f, model_w, model_f) in zip (observedSpectrum.wavelengths, observedSpectrum.flux, s.wavelengths, s.flux):
			chi_sq+= (data_f - model_f)**2
		if chi_sq<lowest_chi:
			lowest_chi = chi_sq
			bestFit = s.objectName
			bestIndex = index
		xValues = s.wavelengths
		yValues = s.flux
		residuals = [of - mf for (of, mf) in zip(observedSpectrum.flux, s.flux)]
		offset = 0.015
		yValues = [y + offset for y in yValues]
		yObserved = [y + offset for y in observedSpectrum.flux]
		matplotlib.pyplot.plot(xValues, yValues,  drawstyle = 'steps', color = 'b')
		matplotlib.pyplot.plot(observedSpectrum.wavelengths, yObserved, color='k', drawstyle = 'steps')
		matplotlib.pyplot.plot(xValues, residuals, color='r', drawstyle = 'steps')
		matplotlib.pyplot.plot([min(xValues), max(xValues)], [0, 0], color='k', linestyle='dashed')
		labelX = 6000
		labelY = observedSpectrum.getNearestFlux(labelX) + 0.01
		print labelX, labelY, s.objectName, chi_sq
		matplotlib.pyplot.text(labelX, labelY, s.objectName + "[%f]"%chi_sq, fontsize=15)

		fig = matplotlib.pyplot.gcf()
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block = False)
		fig.savefig('fit_%s.eps'%s.objectName,dpi=100, format='eps')
		fig.savefig('fit_%s.png'%s.objectName,dpi=200, format='png')

	print "Best fit %s %d.  [%f]"%(bestFit, bestIndex, lowest_chi)

	# Now plot the best fit again...
	matplotlib.pyplot.figure(figsize=(14, 6))
	model = spectra[bestIndex]
	# subSpectrum = copy.deepcopy(model)
	# newLength = subSpectrum.trimWavelengthRange(7000, 7500)
	# area = subSpectrum.integrate()
	# model.divide(area)
	# model.divide(1/tioArea)
	offset = 0.5
	offsetModel = [f + offset for f in model.flux]
	residuals = [f - m for f, m in zip(observedSpectrum.flux, model.flux)]
	matplotlib.pyplot.plot(model.wavelengths, offsetModel,  drawstyle = 'steps', color = 'b')
	matplotlib.pyplot.plot(observedSpectrum.wavelengths, observedSpectrum.flux, color='k', drawstyle = 'steps')
	matplotlib.pyplot.plot(observedSpectrum.wavelengths, residuals, color='r', drawstyle = 'steps')
	matplotlib.pyplot.plot([min(observedSpectrum.wavelengths), max(observedSpectrum.wavelengths)], [0, 0], color='k', linestyle='dashed')
	matplotlib.pyplot.yticks(fontsize = 18)
	matplotlib.pyplot.ylabel(r"F$_{\nu}$ (mJy)", size = 18)
	matplotlib.pyplot.xticks(fontsize = 18)
	matplotlib.pyplot.xlabel(r"Wavelength $(\AA)$", size = 18)
	
	fig = matplotlib.pyplot.gcf()
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block = True)
	fig.savefig('modelfit.eps',dpi=100, format='eps')
	fig.savefig('modelfit.png',dpi=200, format='png')
	
	
	mDwarfTemplate = spectra[bestIndex]
	mDwarfTemplate.writeCSV("mDwarf.csv")


