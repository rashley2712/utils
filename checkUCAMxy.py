#!/usr/bin/env python
import argparse, sys
import numpy
import time
import matplotlib.pyplot

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Tails a UCAM generated log file and keeps track of the X-Y positions of an aperture..')
	parser.add_argument('logfile', type=str, help='ULTRACAM log file')
	parser.add_argument('-a', '--aperture', type=int, default = 1, help='Aperture number. Defaults to 1 if not specified.')
	parser.add_argument('-c','--ccd', type=int, default = 1, help='CCD number. Defaults to 1 if not specified.')
	parser.add_argument('-s','--start', type=int, default = 10, help='Start frame for calculating target position. Defaults to 10.')
	parser.add_argument('-e','--end', type=int, default = 50, help='End frame for calculating target position. Defaults to 50.')
	parser.add_argument('-t','--sleep', type=int, default = 50, help='Sleep time to wait between checking the log file for new data. Defaults to 5 seconds.')
	parser.add_argument('--scale', type=float, default = 0.3, help='Field scale in arcseconds per pixel. Defaults to 0.3 (WHT).')
	 
	arg = parser.parse_args()
	print arg
	
	matplotlib.pyplot.figure(figsize=(12, 8))
	matplotlib.pyplot.xlabel("Frame number", size = 14)
	matplotlib.pyplot.subplot(2, 1, 1)
	matplotlib.pyplot.ylabel("X (pixels)", size = 14)
	matplotlib.pyplot.subplot(2, 1, 2)
	matplotlib.pyplot.ylabel("Y (pixels)", size = 14)
	
	desiredAperture = arg.aperture
	desiredCCD = arg.ccd
	startFrame = arg.start
	endFrame = arg.end
	ccdIndex = 4
	fieldScale = arg.scale
	sleep = arg.sleep	
	
	stop = False
	while not stop:
		inputFile = open(arg.logfile,'r')
		
		
		
		xIndex = 14*(desiredAperture-1) + 8
		yIndex = 14*(desiredAperture-1) + 9
		xValues = []
		yValues = []
		for line in inputFile:
			if line[0] != '#':
				params = line.split()
				CCD = int(params[ccdIndex])
				if CCD == desiredCCD:
					xPosition = float(params[xIndex])
					yPosition = float(params[yIndex])
					# print "%d (%f, %f)"%(CCD, xPosition, yPosition)
					xValues.append(xPosition)
					yValues.append(yPosition)
		
		meanX = numpy.mean(xValues)
		meanY = numpy.mean(yValues)
		sigmaX = numpy.std(xValues)
		sigmaY = numpy.std(yValues)
		print "Total frames read:", len(xValues)
		print "Mean position (%f, %f)"%(meanX, meanY)	
		print "Standard dev  (%f, %f)"%(sigmaX, sigmaY)	
		print "Last position was: (%f, %f)"%(xValues[-1], yValues[-1])
		if len(xValues) > endFrame-1:
			staticX = numpy.mean(xValues[startFrame-1:endFrame-1])
			staticY = numpy.mean(yValues[startFrame-1:endFrame-1])
			print "Target position based on frames %d to %d is (%f, %f)."%(startFrame, endFrame, staticX, staticY)
			deltaX = staticX - xValues[-1]
			deltaY = staticY - yValues[-1]
			deltaXField = deltaX * fieldScale
			deltaYField = deltaY * fieldScale
			print "Drift from target is (%f, %f)  in frame number %d."%(deltaX, deltaY, len(xValues)+1)
			print "Correction to be applied is (%f, %f) arcseconds."%(deltaXField, deltaYField)
		else:
			print "Need more frames to accumulate first, starting at %d"%endFrame
			
		print "------------------------------------------------------------------------------------"

		frames = range(1, len(xValues)+1)
		matplotlib.pyplot.subplot(2, 1, 1)
		matplotlib.pyplot.plot( [frames[0], frames[-1]], [meanX+sigmaX, meanX+sigmaX], color='r', linestyle='dashed')
		matplotlib.pyplot.plot( [frames[0], frames[-1]], [meanX, meanX], color='r', linestyle='-')
		matplotlib.pyplot.plot( [frames[0], frames[-1]], [meanX-sigmaX, meanX-sigmaX], color='r', linestyle='dashed')
		matplotlib.pyplot.plot(frames, xValues, color = 'r')
		
		matplotlib.pyplot.subplot(2, 1, 2)
		matplotlib.pyplot.plot( [frames[0], frames[-1]], [meanY+sigmaY, meanY+sigmaY], color='g', linestyle='dashed')
		matplotlib.pyplot.plot( [frames[0], frames[-1]], [meanY, meanY], color='g', linestyle='-')
		matplotlib.pyplot.plot( [frames[0], frames[-1]], [meanY-sigmaY, meanY-sigmaY], color='g', linestyle='dashed')
		matplotlib.pyplot.plot(frames, yValues, color = 'g')
		
		matplotlib.pyplot.show(block = False)
		
		time.sleep(5)
		
	sys.exit()
