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
			
	def getExtremes(self):
		for index, t in enumerate(self.targets):
			earliest = numpy.min(t['HJD'])
			latest = numpy.max(t['HJD'])
			if index==0:
				minTime = earliest
				maxTime = latest
			if earliest<minTime:
				minTime = earliest
			if latest>maxTime:
				maxTime = latest
		return minTime, maxTime
			
		
	def getNumber(self):
		return len(self.targets)	
		
	def getObs(self, index):
		return self.targets[index]
		
	def getObsByName(self, name):
		for t in self.targets:
			if t['name'] == name: return t
	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Attempt to calculate the probability of detecting RV variability based on the times of observations and a simulated period distribution. ')
	#parser.add_argument('filename', type=str, nargs='+', help='Filename of CSV file..')	
	parser.add_argument('--device', type=str, default='/xs', help = "Dump to the following device. Default is '/xs'.")
	parser.add_argument('--input', type=str, help='File containing a sample of periods to use as a period distribution.')
	parser.add_argument('-p', '--probability', type=str, help='File containing a probability distribution in log P.')
	
	parser.add_argument('-n', type=int, default=100, help='Number of random periods to generate. Default is 100.')
	parser.add_argument('-o', '--observations', type=str, help='Observation times')
	parser.add_argument('--fake', type=int, help='Fake [n] observations evenly spaced over baseline of all observations. Specify [n].')
	
	arg = parser.parse_args()

	
	names = []
	samplePeriods = []
	if arg.input is not None:
		inputFile = open(arg.input, 'rt')
		headers = inputFile.readline().strip()
		
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
		flatDistribution = numpy.random.rand(arg.n) * 4 - 1
		samplePeriods = 10**(flatDistribution)
		logSamplePeriods = numpy.log10(samplePeriods)
		
	periods = samplePeriods
	logPeriods = numpy.log10(periods)	
	# periods = numpy.random.rand(1000)

	if arg.probability is not None:
		inputFile = open(arg.probability, 'rt')
		headers = inputFile.readline().strip()
		bins = []
		probabilities = []
		for line in inputFile:
			fields = line.strip().split('\t')
			print fields
			try:
				bin = float(fields[0])
				probability = float(fields[1])
				bins.append(bin)
				probabilities.append(probability)
			except (ValueError, IndexError):
				print "Could not parse a line in the file"
		binWidth = bins[1] - bins[0]	
		probabilities = numpy.array(probabilities) / sum(probabilities)
		print bins
		print probabilities
		print sum(probabilities)
		generalUtils.query_yes_no("Continue?")
	
		samplePeriods = []	
		for i in range(arg.n):
			testPeriod = numpy.random.choice(bins, p = probabilities)	
			delta = numpy.random.rand() * binWidth
			print testPeriod, delta, testPeriod + delta
			period = 10**(testPeriod + delta)
			print "Sample period: ", period
			samplePeriods.append(period)
			
		probabilityFigure = plt.figure()
		plt.step(bins, probabilities, where='post')
		logPeriods = numpy.log10(samplePeriods)	
		weights = numpy.ones_like(logPeriods)/float(len(logPeriods))
		n, newbins, patches = plt.hist(logPeriods, weights=weights, bins = bins, normed = False, alpha=0.5, color = 'g')
		print zip(n, newbins)
		plt.xlabel('$P_{orb}$ [d]')
		plt.ylabel('Probability')
		plt.title('Input period distribution')
		plt.draw()
		plt.show(block = False)
	
	periods = samplePeriods
	logPeriods = numpy.log10(periods)	
	
	generalUtils.query_yes_no("Continue?")
	

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
	
	K2s = []
	for p, i in zip(periods, inclinations):
		K2s.append(massFunction(0.6, 0.2, p, i))
		print "Period: %f [days], inclination: %f [deg], K2: %f km/s"%(p, i/numpy.pi*180, K2s[-1]) 
		
	phases = numpy.random.rand(len(K2s))
	
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
	
	if arg.fake is not None:
		earliestObs, latestObs = observations.getExtremes()
		obsLength = latestObs - earliestObs
		numObs = arg.fake
		fakeObs = [earliestObs + obsLength/30 * o for o in range(30)]
		observations.addTarget('fake')
		for f in fakeObs:
			observations.addDataToTarget('fake', f)
	
	observations.dumpTargets()
	generalUtils.query_yes_no("Continue?")
	
	def rv(K2, period, phase, date):
		return K2 * numpy.sin(2*numpy.pi/period*date + 2*numpy.pi*phase)
	
	rvSTDs = []
	for period, K2, phase in zip(periods, K2s, phases):
		index = int(numpy.random.rand() * observations.getNumber())
		obs = observations.getObs(index)
		if arg.fake is not None:
			obs = observations.getObsByName('fake')
		print "Random values are: %2.2f days and %.f km/s"%(period, K2)
		print "Chosen observation is number %d which was %s."%(index, obs['name'])
		measuredRVs = []
		for date in obs['HJD']:
			measuredRV = rv(K2, period, phase, date)
			print date, measuredRV
			measuredRVs.append(measuredRV)
			
		rvScatter = numpy.std(measuredRVs)
		
		print "RV scatter is ", rvScatter 
		rvSTDs.append(rvScatter)
	
	detected = []
	not_detected = []
	for p, rv in zip(periods, rvSTDs):
		if rv< 30:
			# print "Not detected", p, rv
			not_detected.append(p)
		else:
			# print "Detected", p, rv
			detected.append(p)
	
	percentageNotDetected = float(len(not_detected))/float(len(detected)) * 100.
	print "Detected %d, Not detected %d  - ratio: %f"%(len(detected), len(not_detected), percentageNotDetected)
	
	figure4 = plt.figure()
	not_detectedDist, bins, patches = plt.hist(numpy.log10(not_detected), 14, facecolor='red', alpha=0.75, normed=False, cumulative=False)
	plt.xlabel('$log_{10}(P_{orb})$ [d]')
	plt.ylabel('N')
	plt.title('Non-detected periods (total sample size: %d)'%arg.n)
	plt.grid(True)
	plt.show(block=False)
	
	figure5 = plt.figure()
	detectedDist, bins, patches = plt.hist(numpy.log10(detected), 14, facecolor='green', alpha=0.75, normed=False, cumulative=False)
	plt.xlabel('$log_{10}(P_{orb})$ [d]')
	plt.ylabel('N')
	plt.title('Detected periods (%d of total sample size: %d)'%(len(detected), arg.n))
	plt.grid(True)
	plt.show(block=False)
	
	index = 0
	
	percentages = []
	probability_not_detection = []
	for n, d in zip(not_detectedDist, detectedDist):
		print n, d, bins[index], bins[index+1], n/d*100
		percentages.append(n/(d))
		probability_not_detection.append(n/(n+d))
		index+=1
		
	figure6 = plt.figure()
	width = bins[1] - bins[0]
	plt.bar(bins[:-1], probability_not_detection, width = width)
	plt.xlabel('$log_{10}(P_{orb})$ [d]')
	plt.ylabel('p(logP)')
	plt.title('Probability of missing an RV detection')
	plt.grid(True)
	plt.show(block=False)
	
	generalUtils.query_yes_no("Continue?")
	
	
