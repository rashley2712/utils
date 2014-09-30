#!/usr/bin/env python

import astropy.io.fits
import argparse
import matplotlib.pyplot, numpy, math
import matplotlib.image
import Image, ImageDraw
import ultracamutils
import classes, wcsclasses
import os, subprocess, sys, json

def getArrayFromObjects(objects, propertyName):
	values = []
	for o in objects:
		value = o[propertyName]
		values.append(value)
	return numpy.array(values)

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Plots the object offset from r, g, b')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-d', '--debuglevel', default = 1, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);

	""" First check if a directory is made for the output files in the working director folder. and create one. 
	"""
	runDate, runNumber = ultracamutils.separateRunNameAndDate(arg.runname)
	
	channels = ['r', 'g', 'b']
	channelDescriptions = {'r':'Red', 'g':'Green', 'b':'Blue'}
	
		
	""" Read info about the run to get the starting RA and DEC coordinates
	"""
	debug.write("Getting run info from the file:" + config.RUNINFO, level = 2)
	runInfo = ultracamutils.getRunInfo(config.RUNINFO, arg.runname)
	
	debug.write("Run Info:\n----------------------", level = 2)
	debug.write(runInfo, level = 2)
	debug.write("----------------------", level = 2)

	
	jsonFilename = ultracamutils.addPaths(config.SITE_PATH, arg.runname) + "_objects.json"
	
	JSONfile = open(jsonFilename, "r")

	wholeFileString = JSONfile.read()

	allObjectsJSON = json.loads(wholeFileString)

	objects = []

	for i in allObjectsJSON:
		ob = json.loads(i)
		object = {}
		object['id'] = ob['id']
		colourIDs = ob['colourID']
		object['colourID'] = colourIDs
		object['meanPosition'] = ob['meanPosition']
		objects.append(object)
		print object['id']
		
	figure1 = matplotlib.pyplot.figure(figsize=(12, 12))
	
	for o in objects:
		if o['colourID']['r']!=-1:
			r_x,r_y = o['meanPosition']['r'][0], o['meanPosition']['r'][1]
			print r_x, r_y
			if o['colourID']['g']!=-1:
				g_x, g_y = o['meanPosition']['g'][0], o['meanPosition']['g'][1]
				matplotlib.pyplot.plot(g_x, g_y, 'g.')
				matplotlib.pyplot.plot([r_x, g_x], [r_y, g_y], lw=1, color='black')
			matplotlib.pyplot.plot(r_x, r_y, 'r.')
			
	
	
	figure2 = matplotlib.pyplot.figure(figsize=(12, 12))
	
	for o in objects:
		if o['colourID']['r']!=-1:
			r_x,r_y = o['meanPosition']['r'][0], o['meanPosition']['r'][1]
			print r_x, r_y
			if o['colourID']['b']!=-1:
				g_x, g_y = o['meanPosition']['b'][0], o['meanPosition']['b'][1]
				matplotlib.pyplot.plot(g_x, g_y, 'b.')
				matplotlib.pyplot.plot(r_x, r_y, 'r.')
				matplotlib.pyplot.plot([r_x, g_x], [r_y, g_y], lw=1, color='black')
	
	matplotlib.pyplot.show()
	
	figure1.savefig('test1.eps',dpi=100, format='eps')
	figure2.savefig('test2.eps',dpi=100, format='eps')
	
	"""
	Now compare the offsets from colour to colour
	"""
	
	"""
	gridSize = 1024
	gridSpacing = 20
	arrowScale = 2.0
	greenOffsetx = numpy.zeros((gridSize, gridSize))
	greenOffsety = numpy.zeros((gridSize, gridSize))
	blueOffsetx = numpy.zeros((gridSize, gridSize))
	blueOffsety = numpy.zeros((gridSize, gridSize))
	
	for i in range(0, gridSize, gridSpacing):
		print i
		for j in range(0, gridSize, gridSpacing): 
			greenPixel = (i, j)
			world = wcsSolutions['g'].getWorldSIP(greenPixel)
			redPixel = wcsSolutions['r'].getPixel(world)
			offset = (greenPixel[0] - redPixel[0], greenPixel[1] - redPixel[1])
			greenOffsetx[i][j] = offset[0]
			greenOffsety[i][j] = offset[1]
			matplotlib.pyplot.plot([i, i+offset[0] * arrowScale], [j, j+offset[1] * arrowScale], lw=1, color='green')
			
			bluePixel = (i, j)
			world = wcsSolutions['b'].getWorldSIP(greenPixel)
			redPixel = wcsSolutions['r'].getPixel(world)
			offset = (bluePixel[0] - redPixel[0], bluePixel[1] - redPixel[1])
			blueOffsetx[i][j] = offset[0]
			blueOffsety[i][j] = offset[1]
			matplotlib.pyplot.plot([i, i + offset[0] * arrowScale ], [j, j + offset[1] * arrowScale ], lw=1, color='blue')
			
		
	print greenOffsetx
	print greenOffsety
			
	print blueOffsetx
	print blueOffsety

	matplotlib.pyplot.show()
	"""

	"""X,Y = numpy.meshgrid( numpy.arange(0,gridSize,1), numpy.arange(0, gridSize, 1) )
	U = numpy.cos(X)
	V = numpy.sin(Y)
	
	#QP = matplotlib.pyplot.quiver(X, Y, greenOffsetx, greenOffsety, scale = 1.5)
	#QP = matplotlib.pyplot.quiver(X, Y, U, V, scale=2)
	#matplotlib.pyplot.quiverkey(QP, """