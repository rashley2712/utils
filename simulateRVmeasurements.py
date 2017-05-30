#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import spectrumClasses, timeClasses, generalUtils
import scipy.optimize
import copy
import ppgplot
import scipy.signal as signal

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

def massFunction(m1, m2, period, i):
	Msol =  1.989e30  		# Solar mass (kg)
	G = 6.678e-11     		# Gravitational constant (m^3/kg/s^2)
	m1 = m1 * Msol
	m2 = m2 * Msol
	p = period * 86400.  	# Period in seconds
	return numpy.power( ((m1 * numpy.sin(i))**3 / (m1 + m2)**2) * 2. * numpy.pi * G / p  , 1./3.) /1000.
	
class sampleObservations:
	def __init__(self):
		self.targets = []
	
	def addTarget(self, name):
		target = {}
		target['name'] = name
		target['HJD'] = []
		self.targets.append(target)
		
	def addDataToTarget(self, targetName, HJD):
		for t in self.targets:
			if t['name'] == targetName: 
				t['HJD'].append(HJD)
				
	def dumpTargets(self):
		for t in self.targets:
			print t['name'], len(t['HJD']), 'observations'
			print t['HJD']	
		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Attempt to calculate the probability of detecting RV variability based on the times of observations and a simulated period distribution. ')
	#parser.add_argument('filename', type=str, nargs='+', help='Filename of CSV file..')	
	parser.add_argument('--device', type=str, default='/xs', help = "Dump to the following device. Default is '/xs'.")
	parser.add_argument('--input', type=str, help='File containing a sample of periods to use as a period distribution.')
	parser.add_argument('-n', type=int, default=100, help='Number of random periods to generate. Default is 100.')
	parser.add_argument('-o', '--observations', type=str, help='Observation times')
	
	arg = parser.parse_args()

	
	names = []
	samplePeriods = []
	if arg.input is not None:
		inputFile = open(arg.input, 'rt')
		for line in inputFile:
			values = line.strip().split('\t')
			names.append(values[0])
			samplePeriods.append(float(values[1]))
		inputFile.close()

		samplePeriods = [p/24. for p in samplePeriods]
		logSamplePeriods = numpy.log10(samplePeriods)
		
		# Draw the histogram of the input data
		figure1 = plt.figure()
		n, bins, patches = plt.hist(logSamplePeriods, 14, facecolor='green', alpha=0.75, cumulative=False)
		print n, bins, patches
		print numpy.sum(n), len(samplePeriods)
		plt.xlabel('$P_{orb}$ [d]')
		plt.ylabel('N')
		plt.title('Nebot period distribution')
		plt.grid(True)

		plt.show(block=False)

		print "Number of period bins", len(n)
		# Normalise the histogram to a probability function	
		probabilities = n / numpy.sum(n)
		binWidth = bins[1] - bins[0]
		print len(bins), len(n)
		print "bin width:", binWidth
		samplePeriods = []
		for i in range(arg.n):
			testPeriod = numpy.random.choice(bins[:-1], p = probabilities)	
			delta = numpy.random.rand() * binWidth
			print testPeriod, delta
			period = 10**(testPeriod + delta)
			print "Sample period: ", period
			samplePeriods.append(period)

	else:
		flatDistribution = numpy.random.rand(arg.n) * 2 - 1
		samplePeriods = 10**(flatDistribution)
	
	periods = samplePeriods
	logPeriods = numpy.log10(periods)	
	# periods = numpy.random.rand(1000)

	figure2 = plt.figure()		
	n, bins, patches = plt.hist(logPeriods, 14, facecolor='grey', alpha=0.75, cumulative=False)
	print n, bins, patches
	print numpy.sum(n), len(samplePeriods)
	plt.xlabel('$P_{orb}$ [d]')
	plt.ylabel('N')
	plt.title('Sample period distribution N = %d'%len(logPeriods))
	plt.grid(True)
	
	plt.show(block=False)
	
	figure3 = plt.figure()
	n, bins, patches = plt.hist(logPeriods, 14, facecolor='green', alpha=0.5, normed=True, cumulative=False)
	n, bins, patches = plt.hist(logSamplePeriods, 14, facecolor='grey', alpha=0.5, normed=True, cumulative=False)
	print n, bins, patches
	print numpy.sum(n), len(samplePeriods)
	plt.xlabel('$P_{orb}$ [d]')
	plt.ylabel('p')
	plt.title('Comparison of period distributions')
	plt.grid(True)
	plt.show(block=False)
	
	# Create the random inclinations. Flat distribution in cos i
	cosi = numpy.random.rand(arg.n)
	inclinations = numpy.arccos(cosi);
	
	""" Plot the inclination probability 
	figure4 = plt.figure()
	n, bins, patches = plt.hist(inclinations*180/numpy.pi, 30, facecolor='green', alpha=0.75, normed=True, cumulative=False)
	plt.xlabel('inclination')
	plt.ylabel('p')
	plt.title('inclination probability')
	plt.grid(True)
	plt.show(block=False)
	"""
	
	
	
	
	K2sini = []
	for p,i in zip(periods, inclinations):
		K2sini.append(massFunction(0.6, 0.2, p, i))
		print "Period: %f [days], inclination: %f [deg], K2sini: %f"%(p, i/numpy.pi*180, K2sini[-1]) 
		
	count = 0
	nonDetectedPeriods = []
	for p, k in zip(periods, K2sini):
		if k<30: 
			count+=1
			nonDetectedPeriods.append(p)
		
	print "%f%% less than 30 km/s"%(float(count)/float(len(K2sini))*100.)

	figure4 = plt.figure()
	n, bins, patches = plt.hist(numpy.log10(nonDetectedPeriods), 14, facecolor='green', alpha=0.75, normed=False, cumulative=False)
	plt.xlabel('$log_{10}(P_{orb})$ [d]')
	plt.ylabel('N')
	plt.title('Non-detected periods')
	plt.grid(True)
	plt.show(block=False)
	
	
	# Load the sample times from the sample times file
	sampleTargets = []
	observations = sampleObservations()
	if arg.observations is not None:
		sampleTimesFile = open(arg.observations, 'rt')
		targetNames = []
		for line in sampleTimesFile:
			items = line.strip().split('\t')
			name = str(items[0])
			HJD = float(items[1])
			if name not in targetNames:
				observations.addTarget(name)
				targetNames.append(name)
			observations.addDataToTarget(name, HJD)
		sampleTimesFile.close()
	
	observations.dumpTargets()
	
	generalUtils.query_yes_no("Continue?")
	
	
	
