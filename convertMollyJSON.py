#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import numpy, math
import argparse
import loadingSavingUtils, statsUtils
import trm.dnl.molly
import spectrumClasses, timeClasses

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a series of spectra that were saved from Molly converts them to JSON format.')
	parser.add_argument('mollyfile', type=str, help='Molly file containing the spectra')
	 
	arg = parser.parse_args()
	print arg

	spectra = []
	
	mollyFilename = arg.mollyfile
	mollyFile = trm.dnl.molly.rmolly(mollyFilename)
	
		
	for index, r in enumerate(mollyFile):
		data = r.toHDU()
		wavelengths = r.x.data
		flux = r.y.data
		head = r.head
		if r.has_mask: print "Mask found"
		spectrum = spectrumClasses.spectrumObject()
		npoints = spectrum.setData(wavelengths, flux)
		targetName = spectrum.parseHeaderInfo(head)
		print "Parsed headers of", targetName
		print r.oneLine()
		spectra.append(spectrum)
		
	numSpectra = len(spectra)

	print "%d spectra loaded."%numSpectra

	for s in spectra:
		outname = "%s_%f.json"%(s.objectName, s.HJD)
		print "Writing to %s"%outname
		s.writeToJSON(outname)
		
