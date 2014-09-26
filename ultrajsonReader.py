#!/usr/bin/env python
import sys, subprocess, re, json, itertools, scipy.stats
import numpy, matplotlib.pyplot 
import argparse
import pylab

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Reads the Ultracam ultra.json file and extracts some run data statistics.')
	
	filename = '/storage/astro1/phsaap/ultracam/logs/ultra.json'
	
	JSONfile = open(filename, "r")

	allObjectsJSON = json.load(JSONfile)

	exposeTimes = []
	allRecords = []
	
	for index, object in enumerate(allObjectsJSON):
		exposeTime = float(object['expose'])
		if ((exposeTime<1000) & (exposeTime>=120)): 
			exposeTimes.append(exposeTime)
			allRecords.append(object)
		else: 
			print "Rejected a record where the exposure time was longer than 1000 minutes or less than 15 minutes"
			print index, object


	#exposeTimes = exposeTimes[0:100]
	print exposeTimes
	exposeTimes = numpy.array(exposeTimes)
	print 'Number of runs:', len(exposeTimes)
	print 'Mean length of run', numpy.mean(exposeTimes)
	print 'Mode length of run', scipy.stats.mode(exposeTimes)
	print "Longest run:", max(exposeTimes), exposeTimes.argmax()
	print allRecords[exposeTimes.argmax()]
	
	binwidth = 10
	bottomofbins = 120
	topofbins = 600
	numbins = 1 + ( (topofbins - bottomofbins) / binwidth)
	bins = numpy.linspace(bottomofbins, topofbins, numbins)
	
	print 'Bins (input)', bins
	
	hist = numpy.histogram(exposeTimes, bins)
	
	print hist
	
	print bins
	
	n, bins, patches = matplotlib.pyplot.hist(exposeTimes, bins, histtype='bar')
	
	fig = matplotlib.pyplot.gcf()
	
	
	matplotlib.pyplot.xlabel('Run length (minutes)')
	matplotlib.pyplot.ylabel('Number of runs')
	
	matplotlib.pyplot.axis([bottomofbins, topofbins, min(n), max(n)])
	
	matplotlib.pyplot.show()
	
	DPI = fig.get_dpi()
	print "DPI:", DPI
	DefaultSize = fig.get_size_inches()
	print "Default size in Inches", DefaultSize
	print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

	
	fig.savefig('hist.eps',dpi=100, format='eps')

