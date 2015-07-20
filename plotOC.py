#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse, sys
import astropy.io.fits
import astropy.stats
import loadingSavingUtils, statsUtils
import scipy.optimize
import time
import random
import timeClasses
import csv
import scipy.optimize


def func(x, a1, a2):
	y = a1 * x + a2
	return y	

def quad(x, a1, a2, a3):
	y = a1 * x * x + a2 *x + a3
	return y

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Plot an O-C.')
	parser.add_argument('filename', type=str, help='CSV file of eclipse times.')
	parser.add_argument('-e', '--ephemerisfilename', type=str, help='Ephemeris .dat file.')
	arg = parser.parse_args()
	print arg
	if arg.filename == None:
		print "Need an ephemeris file."
		sys.exit()

	colourMap = ['b' , 'g', 'r', 'c', 'm', 'y', 'k']
	ephemeris = timeClasses.ephemerisObject()
	
	ephemeris.loadFromFile(arg.ephemerisfilename)
	
	print ephemeris
	
	csvfile = open(arg.filename, 'r')
	reader = csv.reader(csvfile, delimiter=',')
	headings = reader.next()
	columns = []
	for h in headings:
		columnName = h.strip()
		columns.append(columnName)
		
	print "columns are:", columns
	
	loggedCycles = []
	MJDs = []
	errors = []	
	observatories = []
	for line in reader:
		values = [v.strip() for v in line]
		loggedCycle = int(values[0])
		MJD = float(values[1])
		error = float(values[2])
		obs = int(values[3])
		loggedCycles.append(loggedCycle)
		MJDs.append(MJD)
		errors.append(error)
		observatories.append(obs)
	
	print "logged cycles:", loggedCycles
	print "JDs:", MJDs
	
	ocErrors = []
	ocs = []
	cycles = []
	for loggedCycle, MJD, error in zip(loggedCycles, MJDs, errors):
		phase = ephemeris.getPhase(MJD)
		cycle, upper = ephemeris.getOrbits(MJD)
		# print loggedCycle, upper, MJD
		if phase>0.5: phase-=1
		phaseDifference = phase 
		# print phaseDifference
		# if phaseDifference < 0: phaseDifference-=1
		ominusc = phaseDifference * ephemeris.Period
		terror = math.sqrt(error**2 + ephemeris.T0_error**2 + ephemeris.Period_error**2) *86400.
		print upper, MJD, phase, phaseDifference, ominusc, ominusc*86400., error, terror
		ocs.append(ominusc*86400.)
		cycles.append(cycle)
		ocErrors.append(terror)
		
	
	# Map some colours
	oColours = []
	for o in observatories:
		index = o % len(colourMap)
		colour = colourMap[index]
		oColours.append(colour)
	# print oColours
	
	# Break it down by observatory
	allData = []
	for index, o in enumerate(observatories):
		addToExisting = False
		j=0
		for j, a in enumerate(allData):
			if a['obs']== o: addToExisting = True
		
		if not addToExisting:
			data = {}
			data['obs'] = o
			data['colour'] =  colourMap[o % len(colourMap)]
			data['HJDs'] = []
			data['Cycles'] = []
			data['OCs'] = []
			data['OCerrors'] = []
			data['HJDs'].append(MJDs[index])
			data['Cycles'].append(cycles[index])
			data['OCs'].append(ocs[index])
			data['OCerrors'].append(ocErrors[index])
			allData.append(data)
		else:
			allData[j]['HJDs'].append(MJDs[index])
			allData[j]['Cycles'].append(cycles[index])
			allData[j]['OCs'].append(ocs[index])
			allData[j]['OCerrors'].append(ocErrors[index])
			
	
	""" Fit a straight line to TNT points """
	# Retriece just the TNT data...
	tntData = {}
	for d in allData:
		if d['obs'] == 8: tntData = d;
	tntCycles = tntData['Cycles']
	tntOCs = tntData['OCs']
	tntOCErrors = tntData['OCerrors']
	print "TNT Cycles:", tntCycles, len(tntCycles)
	print "TNT OCs:", tntOCs, len(tntOCs)
	print "TNT OC errors:", tntOCErrors, len(tntOCErrors)
	m = 0.
	c=-20.
	guess = numpy.array([m, c])
	x_values = numpy.array(tntCycles)
	y_values = numpy.array(tntOCs)
	y_errors = numpy.array(tntOCErrors)
	print "initial guess", guess
	result = scipy.optimize.curve_fit(func, x_values, y_values, guess, y_errors)
	parameters = result[0]
	m_fit = parameters[0]
	c_fit = parameters[1]
	print "Fit m:%4.8f  c:%4.8f"%(m_fit, c_fit)
	newPeriod = ephemeris.Period + m_fit/86400.
	print "Original period %1.12f  New period %1.12f"%(ephemeris.Period, newPeriod)
	
	# Now fit a straight line to all of the data
	x_values = numpy.array(cycles)
	y_values = numpy.array(ocs)
	y_errors = numpy.array(ocErrors)
	m = 0.
	c=-20.
	guess = numpy.array([m, c])
	print "Now fitting on all of the data points"
	print "initial guess", guess
	result, covariance = scipy.optimize.curve_fit(func, x_values, y_values, guess, y_errors)
	parameters = result
	errors = numpy.sqrt(numpy.diag(covariance))
	print "Covariance:", covariance
	m_err = errors[0]
	m_fit = parameters[0]
	c_fit = parameters[1]
	print "Fit m:%4.8f  c:%4.8f"%(m_fit, c_fit)
	newPeriod = ephemeris.Period + m_fit/86400.
	newError = m_err / 86400.
	newT0 = ephemeris.T0 + c_fit/86400.
	newT0Error = errors[1]/86400.
	print "Original period %1.12f[%e]  New period %1.12f[%e]"%(ephemeris.Period, ephemeris.Period_error, newPeriod, newError)
	print "Original T0 %7.8f[%e]  New T0: %7.8f[%e]"%(ephemeris.T0, ephemeris.T0_error, newT0, newT0Error)
	
		
	""" Try a quadratic. Just for a laugh. """
	x_values = numpy.array(cycles)
	y_values = numpy.array(ocs)
	y_errors = numpy.array(ocErrors)
	# print x_values, y_values, y_errors
	a1 = 0.0
	a2 = 0.0
	a3 = 0.0 
	guess = numpy.array([a1, a2, a3])
	print "initial guess", guess
	result = scipy.optimize.curve_fit(quad, x_values, y_values, guess, y_errors)
	parameters = result[0]
	print "Quadratic result: ", parameters
	a1 = parameters[0]
	a2 = parameters[1]
	a3 = parameters[2]
		
	matplotlib.pyplot.figure(figsize=(16, 8))
	matplotlib.pyplot.xlabel("Cycles", size=22)
	matplotlib.pyplot.ylabel("O-C (seconds)", size=22)
	matplotlib.pyplot.tick_params(axis='both', which='major', labelsize=22)
	#matplotlib.pyplot.scatter(cycles, ocs, c=oColours)
	
	for a in allData:
		cycles = a['Cycles']
		ocs = a['OCs']
		ocErrors = a['OCerrors']
		colour = a['colour']
		matplotlib.pyplot.errorbar(cycles, ocs, color = colour, yerr=ocErrors, fmt = '.', ecolor=colour, capsize=0)
	print cycles
	axes = matplotlib.pyplot.gca()
	xmin, xmax = axes.get_xlim()
	matplotlib.pyplot.plot( [xmin, xmax], [0, 0], color='k', linestyle='dashed')
	#matplotlib.pyplot.plot( [xmin, xmax], [xmin * m_fit + c_fit, xmax * m_fit + c_fit], color='b', linestyle='dashed')
	
	# Draw some special arrows to show our new data points
	for xpos in cycles:
		#matplotlib.pyplot.plot( [ xpos, xpos], [10, 30], color = 'g', linestyle='solid')
		matplotlib.pyplot.arrow( xpos, 30, 0, -20, fc="g", ec="g", head_width=250, head_length=3 )
		
	xPlots = numpy.arange(xmin, xmax, 100)
	yPlots = quad(xPlots, a1, a2, a3)
	# matplotlib.pyplot.plot( xPlots, yPlots, color='g', linestyle='dotted')
	
	matplotlib.pyplot.xlim(xmin, xmax)
	
	fig = matplotlib.pyplot.gcf()
	
	matplotlib.pyplot.show()
	fig.savefig('css081231_oc.eps',dpi=100, format='eps')
	fig.savefig('css081231_oc.png',dpi=100, format='png')
	sys.exit()
	
	print "======================================================="
	print "Red"
	MJD = 57163.17148243
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)
	
	
	
	print "======================================================="
	print "Green"
	MJD = 57163.17148299
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)

	print "======================================================="
	print "Blue"
	MJD = 57163.17148816
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)


	print "======================================================="
	print "LT Rise"
	MJD = 56944.9054808
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)


	print "======================================================="
	print "trm BMJD"
	MJD = 57163.17069173
	
	print "Phase of %f is %f"%(MJD, ephemeris.getPhase(MJD)) 
	print "N Orbits of %f is %d"%(MJD, ephemeris.getOrbits(MJD)) 
	
	phase = ephemeris.getPhase(MJD)
	phaseDifference = 1 - phase
	
	ominusc = phaseDifference * ephemeris.Period
	
	print "O-C is %f days, or %f seconds."%(ominusc, ominusc*86400.)
