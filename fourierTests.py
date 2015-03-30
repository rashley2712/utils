#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import os
import ultracamutils
import matplotlib.pyplot
import argparse
import numpy, math
import classes
import ultraspecClasses
#import rashley_utils as utils
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import ultracam_shift
import time, datetime
import json
from scipy import ndimage
from scipy import fftpack
import Image
import ucamObjectClass
from photutils import datasets
from photutils import daofind
from photutils import aperture_photometry, CircularAperture, psf_photometry, GaussianPSF
import astropy.table, astropy.io
from astropy.stats import median_absolute_deviation as mad


if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Reads a FITS image and plays with Fourier.')
	parser.add_argument('filename', type=str, help='FITS file name.')
	parser.add_argument('--testimage', action='store_true', help='Use the test image.')
	
	arg = parser.parse_args()

	print "About to load the file:", arg.filename
	
	hdulist = astropy.io.fits.open(arg.filename)
	
	print hdulist.info()
	
	for headerItem in hdulist[0].header:
		print headerItem, hdulist[0].header[headerItem]
	
	imageData = hdulist[0].data
	hdulist.close()
	
	if arg.testimage:
		imageData = matplotlib.pyplot.imread("/Users/rashley/code/utils/lena.png")
		print imageData
	
	maximumValue = numpy.max(imageData)
	minimumValue = numpy.min(imageData)
	fullRange = maximumValue - minimumValue
	imageData = imageData/ fullRange * 255.0
	print imageData
	
	originalFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Original image")
	boostedImageData = ultracamutils.percentiles(imageData, 10, 99.8)	
	
	imagePlot = matplotlib.pyplot.imshow(boostedImageData, cmap='gray')
	
	#matplotlib.pyplot.gca().invert_yaxis()			
	matplotlib.pyplot.show(block=False)
	
	freqData = fftpack.fft2(imageData)
	
	transformFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("FFT image")
	
	magnitudes = numpy.absolute(freqData)
	"""print freqData
	print magnitudes
	maximumValue = numpy.max(magnitudes)
	minimumValue = numpy.min(magnitudes)
	fullRange = maximumValue - minimumValue
	magnitudes = magnitudes/ fullRange * 255.0
	print magnitudes
	"""
	magnitudes = ultracamutils.percentiles(magnitudes, 10, 99)
	frequencyPlot = matplotlib.pyplot.imshow(magnitudes, cmap='gray')
	
	matplotlib.pyplot.gca().invert_yaxis()			
	matplotlib.pyplot.show(block=False)
	
	lowPassMask = numpy.ma.zeros(numpy.shape(freqData))
	
	highPassMask = numpy.ma.ones(numpy.shape(freqData))
	
	frequencyCap = 150
	(xmax, ymax) = numpy.shape(lowPassMask)
	lowPassMask[0:frequencyCap, 0:frequencyCap] = 1
	lowPassMask[0:frequencyCap, ymax-frequencyCap:ymax] = 1
	lowPassMask[xmax-frequencyCap:xmax, ymax-frequencyCap:ymax] = 1
	lowPassMask[xmax-frequencyCap:xmax, 0:frequencyCap] = 1
	
	highPassMask[0:frequencyCap, 0:frequencyCap] = 0
	highPassMask[0:frequencyCap, ymax-frequencyCap:ymax] = 0
	highPassMask[xmax-frequencyCap:xmax, ymax-frequencyCap:ymax] = 0
	highPassMask[xmax-frequencyCap:xmax, 0:frequencyCap] = 0
	
	
	filteredFrequency = freqData * highPassMask
	
	
	maskedFFTFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Masked FFT")
	#matplotlib.pyplot.gca().invert_yaxis()			
	maskedFFTPlot = matplotlib.pyplot.imshow(ultracamutils.percentiles(numpy.absolute(filteredFrequency), 10, 99), cmap='gray')
	matplotlib.pyplot.show(block=False)
	
	filteredImage = fftpack.ifft2(filteredFrequency)
	
	
	finalFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Filtered image")
	boostedImageData = ultracamutils.percentiles(numpy.abs(filteredImage), 10, 99.8)	
	
	imagePlot = matplotlib.pyplot.imshow(boostedImageData, cmap='gray')
	
	#matplotlib.pyplot.gca().invert_yaxis()			
	matplotlib.pyplot.show(block=False)
	
	imageData = imageData - numpy.abs(filteredImage)
	
	finalFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Final image")
	boostedImageData = ultracamutils.percentiles(imageData, 10, 99.8)	
	
	imagePlot = matplotlib.pyplot.imshow(boostedImageData, cmap='gray')
	
	#matplotlib.pyplot.gca().invert_yaxis()			
	matplotlib.pyplot.show(block=True)
	
	
	sys.exit()
	
	
	
	# Reconstruct the full frame from the windows	
	stackedFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Initial 10 frame stacked image")
	boostedFullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		boostedImage = ultracamutils.percentiles(w.stackedData, 10, 99.8)
		image = w.stackedData
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		boostedFullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + image
		
		bkg_sigma = 1.48 * mad(image)
		print "bkg_sigma", bkg_sigma   
		sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma) 
		print "Num sources in this window:", len(sources)  
		w.setSourcesAvoidBorders(sources)	
		
		
	# Get the source list from this image
	# Combine the sources in all of the windows
	allSources = []
	for index, w in enumerate(allWindows):
		xll = w.xll/w.xbin - xmin
		yll = w.yll/w.ybin - ymin
		sources = w.getSources()
		positions = zip(sources['xcentroid'], sources['ycentroid'], sources['flux'])
		new_positions = [(x + xll, y + yll, flux) for (x, y, flux) in positions]
		allSources+=new_positions
		
		
	allSources = sorted(allSources, key=lambda object: object[2], reverse = True)
	numSources = len(allSources)
	maxSources = int(round((numSources)*0.6))
	topSources = allSources[0:maxSources]
	print "Number of sources: %d, number of top sources: %d"%(numSources, maxSources)
	# Display the image on the user's screen
	image = matplotlib.pyplot.imshow(boostedFullFrame, cmap='gray_r')
	for s in allSources:
		x, y = s[0], s[1]
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 10, color='green', fill=False, linewidth=1.0))
	for s in topSources:
		x, y = s[0], s[1]
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 10, color='blue', fill=False, linewidth=1.0))
	
	matplotlib.pyplot.gca().invert_yaxis()			
	matplotlib.pyplot.show(block=False)

	masterApertureList = [ (x, y) for (x, y, flux) in topSources]
		
	""" End of the prework """

	rdat.set(1)		# Reset back to the first frame
	
	frameRange = maximumFrames - startFrame + 1
	
	if arg.numframes!=None:
		requestedNumFrames = arg.numframes
		if requestedNumFrames<(frameRange):
			frameRange = requestedNumFrames
	
	startTime = datetime.datetime.now()
	timeLeftString = "??:??"
	""" Run through all the frames in the .dat file.
	"""
	if arg.preview:
		matplotlib.pyplot.figure(figsize=(8, 8))
		matplotlib.pyplot.ion()
		fig = matplotlib.pyplot.gcf()
		matplotlib.pyplot.title("Frame image")
		if arg.stack:
			matplotlib.pyplot.title("Stacked image")
		zoomedImage = matplotlib.pyplot.figure(figsize=(5, 5))
		
			
	fullFrame = numpy.zeros((1057, 1040))
			
	allWindows = []

	ccdFrame = rdat()
	#ccdFrame.rback()
	window = ccdFrame[0]
	
	for windowIndex, w in enumerate(window):
		# Set up some info about the window sizes and extents
		window = ultraspecClasses.window()
		window.setExtents(w.llx, w.lly, w.nx, w.ny)
		window.setBinning(w.xbin, w.ybin)
		
		image = w._data
		#image -= numpy.median(image)
		if (arg.usefirstframe):
			window.setData(image)
			bkg_sigma = 1.48 * mad(image)
			print "bkg_sigma", bkg_sigma   
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			window.setSourcesAvoidBorders(sources)	
		else: 
			window.setBlankData(image)
		
		
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin
	sourceMap = ultraspecClasses.sourceMap((fullFrameysize, fullFramexsize))
			
	debug.write("Building a map of sources in order to define the apertures...", level = 2)
	for frameIndex in range(2, frameRange + 1):
		framesToGo = frameRange - frameIndex
		currentTime = datetime.datetime.now()
		trueFrameNumber = startFrame + frameIndex - 1
		completionPercent = (float(frameIndex) / float(frameRange) * 100.)
		timePassed = ultracamutils.timedeltaTotalSeconds(currentTime - startTime)
		totalTime = timePassed * 100. / completionPercent
		etaTime = startTime + datetime.timedelta(seconds = totalTime)
		timeLeft = etaTime - currentTime
		(hours, mins, secs) = ultracamutils.timedeltaHoursMinsSeconds(timeLeft)
		timeLeftString = str(hours).zfill(2) + ":" + str(mins).zfill(2) + ":" + str(secs).zfill(2)
		
		ccdFrame = rdat()
		
		statusString = "\r%s Frame: [%d/%d]"%(timeLeftString, trueFrameNumber, frameRange)
		sys.stdout.write(statusString)
		sys.stdout.flush()
		
		windows = ccdFrame[0]
		
		for windowIndex, w in enumerate(windows):
			image = w._data
			#image -= numpy.median(image)
			#allWindows[windowIndex].addData(image)
			bkg_sigma = 1.48 * mad(image)
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			#print sources
			allWindows[windowIndex].setSourcesAvoidBorders(sources)	
			
		# Combine the sources in all of the windows
		allSources = []
		for index, w in enumerate(allWindows):
			xll = w.xll/w.xbin - xmin
			yll = w.yll/w.ybin - ymin
			sources = w.getSources()
			positions = zip(sources['xcentroid'], sources['ycentroid'], sources['flux'])
			new_positions = [(x + xll, y + yll, flux) for (x, y, flux) in positions]
			allSources+=new_positions

		allSources = sorted(allSources, key=lambda object: object[2], reverse = True)
		#print allSources
		tempSources = [ (x, y) for (x, y, flux) in allSources]
		#print tempSources
		allSources = tempSources
		sourceMap.updateMap(allSources)
		
		if (applyShift):
		
			oldCatalog = numpy.array(masterApertureList)
			newCatalog = numpy.array(allSources)

			psize  = 0.1
			fwhm   = 4.
			dmax   = 10.
			mmax   = 10.

			(gaussImage, xp, yp, xr, yr) = ultracam_shift.vimage(oldCatalog, newCatalog, dmax, psize, fwhm)
			# (nmatch, inds) = ultracam_shift.match(oldCatalog, newCatalog, xp, yp, mmax)
			debug.write("Calculated offset: (%2.2f, %2.2f)"%(xr, yr), level = 2)

		for windowIndex, w in enumerate(windows):
			image = w._data
			allWindows[windowIndex].setData(image)
			if (applyShift): image = ndimage.interpolation.shift(image, (-1.0*yr, -1.0*xr), order = 3 )
			allWindows[windowIndex].addToStack(image)
			
			
		if arg.preview: 
			fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
			for w in allWindows:
				if (arg.stack):
					boostedImage = ultracamutils.percentiles(w.stackedData, 20, 99)
				else:
					boostedImage = ultracamutils.percentiles(w.data, 20, 99)
				xll = w.xll/w.xbin - xmin
				xsize = w.nx
				yll = w.yll/w.ybin - ymin
				ysize = w.ny
				fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
					
			matplotlib.pyplot.figure(fig.number)
			matplotlib.pyplot.imshow(fullFrame, cmap='gray_r')
			
			for s in allSources:
				(x, y) = s
				matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 15, color='green', fill=False, linewidth=1.0))
			
			
			matplotlib.pyplot.title("Frame image [%d/%d]"%(trueFrameNumber, frameRange))
			if arg.stack:
				matplotlib.pyplot.title("Stacked image [%d/%d]"%(trueFrameNumber, frameRange))
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show()
			matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'
			
			# Now also draw the zoomed in region around the first aperture
			(x, y) = masterApertureList[0]
			croppedFrame = fullFrame[y-10:y+10, x-9:x+10]
			
			matplotlib.pyplot.figure(zoomedImage.number)
			matplotlib.pyplot.imshow(croppedFrame, cmap='gray_r', interpolation = 'none')
			matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((10,10), 1, color='green', fill=False, linewidth=1.0))
			if (applyShift): matplotlib.pyplot.plot([10, 10+xr], [ 10, 10+yr], lw=1, color='green')
			matplotlib.pyplot.title("Zoom on aperture number 1: Frame [%d/%d]"%(trueFrameNumber, frameRange))
			matplotlib.pyplot.xlim([0, 20])
			matplotlib.pyplot.ylim([0, 20])
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show()
			
			matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'
			
		if arg.sleep!=0:
			time.sleep(arg.sleep)
	sys.stdout.write("\rProcessed %d frames      \n"%frameRange)
	sys.stdout.flush()
	
	
	# Draw the final stacked frame and it's sources
	matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.ion()
	finalFigure = matplotlib.pyplot.gcf()
	fullFrameImage = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		boostedImageData = ultracamutils.percentiles(w.stackedData, 20, 99)
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrameImage[yll:yll+ysize, xll:xll+xsize] = fullFrameImage[yll:yll+ysize, xll:xll+xsize] + boostedImageData
	
	matplotlib.pyplot.title("Final stacked image")
	matplotlib.pyplot.imshow(fullFrameImage, cmap='gray_r')
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block = True)
	matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'
	
	
	
	# Get the source map
	smoothedSourceMap = sourceMap.getSmoothMap()
	
	# Now use this source map to generate a set of apertures
	bkg_sigma = 1.48 * mad(smoothedSourceMap)
	print "sourceMap median:", numpy.median(smoothedSourceMap)
	print "sourceMap mean:", numpy.mean(smoothedSourceMap)
	print "sourceMap max:", numpy.max(smoothedSourceMap)
	threshold = frameRange/100.
	print "threshold:", threshold
	apertureSources = daofind(smoothedSourceMap, fwhm=4.0, threshold=threshold)   
	
	# Draw the source map 
	sourceMapImage = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Source map")
	matplotlib.pyplot.imshow(smoothedSourceMap, cmap='hot')
	for s in apertureSources:
		x, y = s['xcentroid'], s['ycentroid']
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 10, color='green', fill=False, linewidth=1.0))
	matplotlib.pyplot.gca().invert_yaxis()			
	#matplotlib.pyplot.show(block=False)
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sourcemap.png"
	matplotlib.pyplot.savefig(outputFilename)

	# Draw the source map with no apertures
	sourceMapImage = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Source map")
	matplotlib.pyplot.imshow(smoothedSourceMap, cmap='hot')
	matplotlib.pyplot.gca().invert_yaxis()
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sourcemap_clean.png"
	matplotlib.pyplot.savefig(outputFilename)


	# Write the XYLS FITS file
	if (arg.xyls):
		IDs = []
		x_values = []
		y_values = []
		fluxes = []
		sortedObjects = sorted(apertureSources, key= lambda p:p['flux'], reverse=True)
		
		for num, s in enumerate(sortedObjects):
			IDs.append(num)
			x_values.append(s['xcentroid'])
			y_values.append(s['ycentroid'])
			fluxes.append(s['flux'])
			

		FITSFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sources.xyls"
		debug.write("Writing FITS file: " + FITSFilename, level=2)
		col1 = astropy.io.fits.Column(name='ID', format='I', array=IDs)
		col2 = astropy.io.fits.Column(name='X', format='E', array=x_values)
		col3 = astropy.io.fits.Column(name='Y', format='E', array=y_values)
		col4 = astropy.io.fits.Column(name='FLUX', format='E', array=fluxes)
		cols = astropy.io.fits.ColDefs([col1, col2, col3, col4])	
		tbhdu =astropy.io.fits.new_table(cols)
		
		prihdr = astropy.io.fits.Header()
		prihdr['TARGET'] = runInfo.target
		prihdr['RA'] = runInfo.ra
		prihdr['DEC'] = runInfo.dec
		prihdr['COMMENT'] = "This file created by uspecCreateSourceMap.py from the Ultracam pipeline."
		prihdr['RUNIDENT'] = arg.runname
		prihdr['WIDTH'] = fullFramexsize
		prihdr['HEIGHT'] = fullFrameysize
		
		prihdu = astropy.io.fits.PrimaryHDU(header=prihdr)
		thdulist = astropy.io.fits.HDUList([prihdu, tbhdu])
		thdulist.writeto(FITSFilename, clobber=True)
	
	# Generate the stacked image for writing to disc
	stackedFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Stacked image")
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		boostedImage = ultracamutils.percentiles(w.stackedData, 40, 99.8)
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
	
	#image = matplotlib.pyplot.imshow(fullFrame, cmap='gray_r')
	#matplotlib.pyplot.gca().invert_yaxis()			
	#matplotlib.pyplot.show(block=True)
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + ".png"
	#matplotlib.pyplot.savefig(outputFilename)
	
	# Do the same thing, but use the PIL library
	imgData = numpy.rot90(fullFrame, 3)
	imgSize = numpy.shape(imgData)
	imgLength = imgSize[0] * imgSize[1]
	testData = numpy.reshape(imgData, imgLength, order="F")
	img = Image.new("L", imgSize)
	palette = []
	for i in range(256):
		palette.extend((i, i, i)) # grey scale
		img.putpalette(palette)
		
	img.putdata(testData)
	debug.write("Writing PNG file: " + outputFilename, level = 2) 
	img.save(outputFilename, "PNG", clobber=True)
	
	palette = []
	for i in range(256):
		palette.extend((255-i, 255-i, 255-i)) # inverse grey scale
		img.putpalette(palette)
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_inverted.png"
	debug.write("Writing PNG file: " + outputFilename, level = 2) 
	img.save(outputFilename, "PNG", clobber=True)
	
	# Write out the stacked image as a non-normalised FITS image
	FITSFilename =  ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_stacked.fits"
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		imageData = w.stackedData
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + imageData

	print "Writing FITS image"
	median = numpy.median(fullFrame)
	maximum = numpy.max(fullFrame)
	minimum = numpy.min(fullFrame)
	print "Median of the image:", median
	print "Maximum of the image:", maximum
	print "Minimum of the image:", minimum
	fullFrame = fullFrame - median
	fullFrame = numpy.clip(fullFrame, 0, 1E14)
	bkg_sigma = 1.48 * mad(fullFrame)
	print "bkg_sigma", bkg_sigma   
	sources = daofind(fullFrame, fwhm=4.0, threshold=3*bkg_sigma)   
	
	print "Final sources:"
	print sources
	
	
	ra = runInfo.ra  # Convert RA to degrees
	dec = runInfo.dec
	fieldScaleX = -8.3E-05
	fieldScaleY = 8.3E-05
	
	prihdr = astropy.io.fits.Header()
	prihdr['COMMENT'] = "This file created by the Ultracam pipeline."
	prihdr['TARGET'] = runInfo.target
	prihdr['COMMENT'] = runInfo.comment
	prihdr['EQUINOX'] = 2000
	prihdr['RADECSYS'] = "FK5"
	prihdr['CTYPE1'] = "RA---TAN"
	prihdr['CTYPE2'] = "DEC--TAN"
	prihdr['CRPIX1'] = fullFramexsize/2
	prihdr['CRPIX2'] = fullFrameysize/2
	prihdr['CRVAL1'] = ra
	prihdr['CRVAL2'] = dec
	prihdr['CDELT1'] = fieldScaleX
	prihdr['CDELT2'] = fieldScaleY
	
	hdu = astropy.io.fits.PrimaryHDU(fullFrame, header=prihdr)
	
	hdulist = astropy.io.fits.HDUList([hdu])
	hdulist.writeto(FITSFilename, clobber=True)
	
	# Now write out the aperture data
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_apertures.json"
	outputFile = open(outputFilename, "w")
	apertureList = []
	for i, s in enumerate(apertureSources):
		aperture = {}
		aperture['id'] = i
		aperture['x'] = s['xcentroid']
		aperture['y'] = s['ycentroid']
		aperture['sharpness'] = s['sharpness']
		aperture['roundness1'] = s['roundness1']
		aperture['roundness2'] = s['roundness2']
		aperture['flux'] = s['flux']
		apertureList.append(aperture)
	json.dump(apertureList, outputFile)
	outputFile.close()
		
