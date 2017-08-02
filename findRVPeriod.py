#!/usr/bin/env python
import sys, os
import numpy, math
import argparse
import loadingSavingUtils, generalUtils
import spectrumClasses, timeClasses
import scipy.optimize
import copy
import ppgplot
import gSheets
import matplotlib.pyplot

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def findAliases(n, f, p):
	a = []
	lp = copy.deepcopy(p)
	for index in range(n):
		a.append(f[numpy.argmax(lp)])
		lp[numpy.argmax(lp)] = 0
		print index+1, a[-1]
	return a
	
def sine(x, a0, a1):
    y = a0 * numpy.sin(2.*math.pi*(x + a1))
    return y
    
def sineFixedFreq(x, gamma, amplitude, phase):
    global frequency
    f = frequency
    y = gamma + amplitude * numpy.sin(2. * numpy.pi * (f * x - t0)) 
    return y
    
def sineFixedGammaAmplitude(x, period, phase):
    global gamma, amplitude
    w = 2. * numpy.pi / period
    y = gamma + amplitude * numpy.sin(w * x + phase)
    return y
    
def sineFit(x, gamma, k2, freq, phase):
	y = gamma + k2 * numpy.sin(2.*numpy.pi * (x * freq + phase))
	return y
    
def sineFixedFreqFit(x, gamma, k2, phase):
	global bestFrequency
	freq = bestFrequency 
	y = gamma + k2 * numpy.sin(2.*numpy.pi * (x * freq + phase))
	return y

def sineFreqPhase(x, gamma, amplitude, f, p):
    y = gamma + (amplitude * numpy.sin(2.*math.pi*f*(x + p)))
    return y
    
   
def sinePhase(x, gamma, amplitude, phase):
    y = gamma + amplitude * numpy.sin(2. * numpy.pi * x + phase) 
    return y

if __name__ == "__main__":  
	parser = argparse.ArgumentParser(description='Loads a CSV file containing HJDs and RVs and tries to fit a period and sinusoid to the data.')
	# parser.add_argument('inputfile', type=str, help='Filename of the CSV file containing the RVs.')
	parser.add_argument('--plo', type=float, default= 0.05, help='Period (days) to start the search. Default is 0.01 days.') 
	parser.add_argument('--phi', type=float, default= 100.0, help='Period (days) to stop the search. Default is 10 days.') 
	parser.add_argument('-n', '--nalias', type=int, default=1, help='Alias number for the fit. Default is 1.')
	parser.add_argument('objectname', type=str, help='Object name.')
	arg = parser.parse_args()
	# print arg
	plo = arg.plo
	phi = arg.phi
    
	# Load the fitted wavelength data from the Google Doc
	docInstance = gSheets.gSheetObject()
	docInstance.initCredentials()
	docInstance.setDocID('11fsbzSII1u1-O6qQUB8P0RzvJ8MzC5VHIASsZTYplXc')
	docInstance.setObjectName(arg.objectname)
	docInstance.loadAllReadings()
    
	data = docInstance.readings
	print "Loaded data for %s from the Google Doc."%arg.objectname
	print "%d data points loaded."%len(data)
	dates = []
	velocities = []
	velErrors = []
	good = []
	print "HJD\t\tVelocity (km/s)\tVel error"
	for index, d in enumerate(data):
		good = d['good']
		if good==0: 
			print bcolors.WARNING + "%f\t%f\t%f"%(d['HJD'], d['RV'], d['RV error']) + bcolors.ENDC 
		else: 
			print "%f\t%f\t%f"%(d['HJD'], d['RV'], d['RV error'])
			dates.append(d['HJD'])
			velocities.append(d['RV'])
			velErrors.append(d['RV error'])
            
	minimumFrequency = 1/phi
	maximumFrequency = 1/plo
	
	from astropy.stats import LombScargle
	frequency, power = LombScargle(dates, velocities, velErrors, fit_mean = True).autopower(minimum_frequency = minimumFrequency, maximum_frequency = maximumFrequency, samples_per_peak=15)
	# frequency, power = LombScargle(dates, velocities, velErrors).autopower()
	print len(frequency), "points in the periodogram"
	
	bestFrequency = frequency[numpy.argmax(power)]
	print "Best frequency: %f cycles per day"%bestFrequency
	bestPeriod = 1/bestFrequency
	print "%s Best period: %f days or %f hours"%(arg.objectname, bestPeriod, bestPeriod * 24.)
	
	aliases = findAliases(10, frequency, power)
	
	if arg.nalias!=1:
		bestFrequency = aliases[arg.nalias-1]
		print "Choosing alias %d ... frequency: %f cycles/day."%(arg.nalias, bestFrequency)
	
	generalUtils.setMatplotlibDefaults()
	
	matplotlib.pyplot.plot(frequency, power, linewidth=1.0, color='k')
	for a in aliases:
		matplotlib.pyplot.plot([a, a], [0, 1], linestyle='--')
	
	matplotlib.pyplot.show(block = False)
	
	x_values = numpy.array(dates)
	y_values = numpy.array(velocities)
	y_errors = numpy.array(velErrors)
		
	k2 = (max(y_values) - min(y_values)) / 2.0
	gamma = numpy.mean(y_values)
	phase = 0.0
	guess = [gamma, k2, bestFrequency, phase]
	print "Guess:", guess
	upperBounds = [300, 400, maximumFrequency, 1.0 ]
	lowerBounds = [-300, 0 , minimumFrequency, 0.0 ]
	bounds = (lowerBounds, upperBounds)
	
	results, covariance = scipy.optimize.curve_fit(sineFit, x_values, y_values, p0 = guess, sigma = y_errors, bounds=bounds)
	errors = numpy.sqrt(numpy.diag(covariance))
	print "Result:", results
	print "Errors:", errors
	gammaFit = results[0]
	gammaError = errors[0]
	k2Fit = results[1]
	k2Error = errors[1]
	frequencyFit = results[2]
	frequencyError = errors[2]
	phaseFit = results[3]
	phaseError = errors[3]
	periodFit = 1.0/frequencyFit
	periodError = frequencyError/frequencyFit/frequencyFit
	t0 = - phaseFit / (frequencyFit) 
	print "Result of curve fit: "
	print "\tgamma velocity: \t%s[%s] km/s"%generalUtils.formatValueError(gammaFit, gammaError)
	print "\tk2 amplitude: \t%s[%s] km/s"%generalUtils.formatValueError(k2Fit, k2Error)
	print "\tfrequency: \t%s[%s] /d"%generalUtils.formatValueError(frequencyFit, frequencyError)
	print "\tperiod: \t%s[%s] days"%generalUtils.formatValueError(periodFit, periodError)
	print "\tphase: \t\t\t%s[%s]"%generalUtils.formatValueError(phaseFit, phaseError)
	print "T0:", t0
	
	# Calculate phases
	phases = []
	for d in dates:
		phase = ((d - t0) % periodFit)/periodFit 
		phases.append(phase)
	addphases = []
	for p in phases:
		addphases.append(p + 1.0)
	phases.extend(addphases)
	velocities.extend(velocities)
	velErrors.extend(velErrors)
	
	phasedFoldedLightCurve = matplotlib.pyplot.figure()
	matplotlib.pyplot.errorbar(phases, velocities, color='k', yerr=velErrors, fmt='.', ecolor='0.75', capsize=0)

	curve = numpy.arange(0, 2, 0.01)
	curveFit = gammaFit + k2Fit * numpy.sin(2*math.pi*curve)
	matplotlib.pyplot.plot(curve, curveFit, color='k', linewidth=1.0)
	matplotlib.pyplot.plot([0,2], [gammaFit, gammaFit], linewidth=1.0, linestyle='--', color='k', alpha=0.75)
	matplotlib.pyplot.show(block = False)
	
	generalUtils.query_yes_no("Continue?")
	
	
	# Refit gamma, phase and k2 with fixed frequency
	
	k2 = k2Fit
	gamma = gammaFit
	phase = phaseFit
	guess = [gamma, k2, phase]
	print "2nd pass guess:", guess
	upperBounds = [300, 400, 1.0 ]
	lowerBounds = [-300, 0 , 0.0 ]
	bounds = (lowerBounds, upperBounds)
	
	results, covariance = scipy.optimize.curve_fit(sineFixedFreqFit, x_values, y_values, p0 = guess, sigma = y_errors, bounds=bounds)
	errors = numpy.sqrt(numpy.diag(covariance))
	print "Result:", results
	print "Errors:", errors
	gammaFit = results[0]
	gammaError = errors[0]
	k2Fit = results[1]
	k2Error = errors[1]
	phaseFit = results[2]
	phaseError = errors[2]
	t0 = - phaseFit / frequencyFit 
	print "Result of curve fit: "
	print "\tgamma velocity: \t%s[%s] km/s"%generalUtils.formatValueError(gammaFit, gammaError)
	print "\tk2 amplitude: \t%s[%s] km/s"%generalUtils.formatValueError(k2Fit, k2Error)
	print "\t*frequency: \t%s[%s] /d"%generalUtils.formatValueError(frequencyFit, frequencyError)
	print "\t*period: \t%s[%s] days"%generalUtils.formatValueError(periodFit, periodError)
	print "\tphase: \t\t\t%s[%s]"%generalUtils.formatValueError(phaseFit, phaseError)
	print "T0:", t0
	# Calculate phases
	phases = []
	for d in dates:
		phase = ((d-t0) % periodFit)/periodFit 
		phases.append(phase)
	addphases = []
	for p in phases:
		addphases.append(p + 1.0)
	phases.extend(addphases)

	# Now do a stacked plot of pgram and folded RV curve
	
	stackedPlot = matplotlib.pyplot.figure()
	ax1 = matplotlib.pyplot.subplot(2, 1, 1)
	 
	matplotlib.pyplot.errorbar(phases, velocities, color='k', yerr=velErrors, fmt='.', ecolor='0.75', capsize=0)

	curve = numpy.arange(0, 2, 0.01)
	curveFit = gammaFit + k2Fit * numpy.sin(2*math.pi*curve)
	matplotlib.pyplot.plot(curve, curveFit, color='k', linewidth=1.0)
	matplotlib.pyplot.plot([0,2], [gammaFit, gammaFit], linewidth=1.0, linestyle='--', color='k', alpha=0.75)
	
	ax2 = matplotlib.pyplot.subplot(2, 1, 2)
	matplotlib.pyplot.plot(frequency, power, linewidth=1.0, color='k', alpha=0.75)
	zoomFraction = 0.1
	matplotlib.pyplot.plot([bestFrequency*(1-zoomFraction), bestFrequency*(1-zoomFraction)], [0, 1], linestyle='--')
	matplotlib.pyplot.plot([bestFrequency*(1+zoomFraction), bestFrequency*(1+zoomFraction)], [0, 1], linestyle='--')
	
	left, bottom, width, height = [0.75, 0.25, 0.14, 0.2]
	ax3 = stackedPlot.add_axes([left, bottom, width, height])
	zoomedFrequency = []
	zoomedPower = []

	for index, f in enumerate(frequency):
		if (f > bestFrequency*(1-zoomFraction)) and (f < bestFrequency*(1+zoomFraction)):
			zoomedFrequency.append(f)
			zoomedPower.append(power[index])
			
	ax3.plot(zoomedFrequency, zoomedPower, linewidth=1.0, color='k', alpha=0.75)
	ax3.plot([bestFrequency, bestFrequency], [0, 1], linestyle='--')
	matplotlib.pyplot.yticks(visible=False)
	matplotlib.pyplot.xticks(visible=False)
		
	matplotlib.pyplot.show(block = False)
	
	generalUtils.query_yes_no("Continue?")
	
	
	sys.exit()
     
	"""
	
    guess = [gamma, amplitude, phase]
    upperBounds = [100, 200, 1]
    lowerBounds = [-100, 0 , 0]
    bounds = (lowerBounds, upperBounds)
    print "Guess for first fit: ", guess
    results, covariance = scipy.optimize.curve_fit(sinePhase, x_values, y_values, p0 = guess, sigma = y_errors, bounds = bounds)
    errors = numpy.sqrt(numpy.diag(covariance))
        
    gammaFit = results[0]
    gammaError = errors[0]
    amplitudeFit = results[1]
    amplitudeError = errors[1]
    phaseFit = results[2]
    phaseError = errors[2]
    print "Result of curve fit: "
    print "\tgamma velocity: \t%s[%s] km/s"%generalUtils.formatValueError(gammaFit, gammaError)
    print "\tv sin i amplitude: \t%s[%s] km/s"%generalUtils.formatValueError(amplitudeFit, amplitudeError)
    print "\tphase: \t\t\t%s[%s]"%generalUtils.formatValueError(phaseFit, phaseError)
    
    xFit = numpy.arange(0, 2, 0.02)
    yFit = gammaFit + amplitudeFit * numpy.sin(2*numpy.pi*(xFit + phaseFit))
    lc = ppgplot.pgqci()
    ls = ppgplot.pgqls()
        
    ppgplot.pgsci(2)
    ppgplot.pgline(xFit, yFit)
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([0, 2], [gammaFit, gammaFit])
    ppgplot.pgsci(lc)
    ppgplot.pgsls(ls)
    
    gamma = gammaFit
    amplitude = amplitudeFit
    phase = phaseFit
    period = periodDays
    loop = True
    solutionFound = False
    while loop:
        print "Current period = %s%f%s days or %f hours."%(bcolors.BOLD, period, bcolors.ENDC, period*24)
        print "Choice: [A] - Accept current period for final least squares fit."
        print "\t[Z] - Zoom in the periodogram."
        print "\t[R] - Reset the periodogram."
        print "\t[X] - Exit."
        choice = sys.stdin.read(1)
        if choice == 'x': 
            loop = False
            continue
        if choice == 'z':
            print "Choose a period to zoom in by clicking on the periodogram."
            ppgplot.pgslct(pgramPGPlotWindow)
            (x1, y, char) = ppgplot.pgcurs(0, 0)
            print "Zooming in on range %f days to ...."%x1, 
            (x2, y, char) = ppgplot.pgband(4, 0, x1, y)
            range = phi-plo
            print x2, "days."
            period = getPeriodogram(pgramPGPlotWindow, reducedDates, velocities, x1, x2)
            (gamma, amplitude, phase) = fitPhasePlot(phasePGPlotWindow, period, reducedDates, velocities, velErrors)
        if choice == 'r':
            period = getPeriodogram(pgramPGPlotWindow, reducedDates, velocities, arg.plo, arg.phi)
            (gamma, amplitude, phase) = fitPhasePlot(phasePGPlotWindow, period, reducedDates, velocities, velErrors)
        if choice == 'a':
            fullResults = finalFit(period, gamma, amplitude, phase, reducedDates, velocities, velErrors)
            print fullResults
            period, periodError = fullResults[0]
            gamma, gammaError = fullResults[1]
            amplitude, amplitudeError = fullResults[2]
            phase, phaseError = fullResults[3]
            plotData = { 'period': period, 
                         'periodError': periodError,
                         'gamma' : gamma,
                         'gammaError' : gammaError, 
                         'amplitude' : amplitude, 
                         'amplitudeError' : amplitudeError, 
                         'phase' : phase,
                         'phaseError' : phaseError, 
                         'x': reducedDates,
                         'y': velocities, 
                         'yerrors': velErrors,
                         'object': arg.objectname
                        }
            phasePlot(plotData)
            t0 = - 2 * numpy.pi * phase / period 
            t0Error = abs(t0 * numpy.sqrt( (periodError/period)**2 + (phaseError/phase)**2 ))
            t0+= xStart
            print "t0", t0, t0Error
            solutionFound = True
            
        
    if solutionFound:
        logFile = open("findRVlog.csv", 'a')
        logFile.write("%s, %f, %f, %f, %f, %f, %f, %f, %f\n"%(arg.objectname, period, periodError, gamma, gammaError, amplitude, amplitudeError, t0, t0Error))
        logFile.close()        
	
	if outputPS and solutionFound:
		print "Writing PS to ", psFilename
		phasePlot(plotData, device=psFilename + "/ps")
        
    
    ppgplot.pgend()
    sys.exit()
		"""
   
