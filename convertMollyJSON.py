#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import numpy, math
import argparse
import loadingSavingUtils, statsUtils
import trm.molly
import spectrumClasses

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a series of spectra that were saved from Molly converts them to JSON format.')
	parser.add_argument('mollyfile', type=str, help='Molly file containing the spectra')
	parser.add_argument('--suffix', type=str, help='Suffix to add to the end of the filenames.')
	parser.add_argument('-e', type=str, help='[Optional] Ephemeris data file')
	
	 
	arg = parser.parse_args()
	print arg

	if arg.e!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = timeClasses.ephemerisObject()
		ephemeris.loadFromFile(arg.e)
		print ephemeris
	else:
		hasEphemeris = False


	spectra = [] 
	
	mollyFilename = arg.mollyfile
	mollyFile = trm.molly.rmolly(mollyFilename)
		
	for index, r in enumerate(mollyFile):
		wavelengths = []
		flux = []
		fluxErrors = []	
		for f, fe, w in zip(r.f, r.fe, r.wave):
			print w, f, fe
			wavelengths.append(w)
			flux.append(f)
			fluxErrors.append(fe)
 
		head = r.head
		spectrum = spectrumClasses.spectrumObject()
		npoints = spectrum.setData(wavelengths, flux, fluxErrors)
		targetName = spectrum.parseHeaderInfo(head)
		spectrum.wavelengthUnits = "\\A"
		spectrum.fluxLabel = r.label
		spectrum.fluxUnits = r.units
		# spectrum.fluxUnits = "relative counts"
		
		print "Parsed headers of %s for HJD: %f"%(targetName, spectrum.HJD)
		spectra.append(spectrum)
		
	numSpectra = len(spectra)

	print "%d spectra loaded."%numSpectra

	for s in spectra:
		outname = "%s_%f.json"%(s.objectName, s.HJD)
		if arg.suffix!=None:
			outname = "%s_%f_%s.json"%(s.objectName, s.HJD, arg.suffix)
		if hasEphemeris:
			outname = "%s_%f_%s.json"%(s.objectName, ephemeris.getPhase(s.HJD), arg.suffix)
			
		print "Writing to %s"%outname
		s.writeToJSON(outname)
		
