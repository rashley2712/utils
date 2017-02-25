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

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

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
    
    
def sineFreqPhase(x, gamma, amplitude, f, p):
    y = gamma + (amplitude * numpy.sin(2.*math.pi*f*(x + p)))
    return y
    
   
def sinePhase(x, gamma, amplitude, phase):
    y = gamma + amplitude * numpy.sin(2. * numpy.pi * (x + phase)) 
    return y

def getPeriodogram(pgplotHandle, xdata, ydata, plo, phi):
    x = numpy.array(xdata)
    # Subtract the mean from the y-data
    y_mean = numpy.mean(ydata)
    y = numpy.array(ydata - y_mean)
    periods = numpy.linspace(plo, phi, 8000)
    ang_freqs = 2 * numpy.pi / periods
    power = signal.lombscargle(x, y, ang_freqs)
    # normalize the power
    N = len(x)
    power *= 2 / (N * y.std() ** 2)
    print power
    ppgplot.pgslct(pgplotHandle)
    ppgplot.pgeras()
    ppgplot.pgenv(min(periods), max(periods), 0, max(power), 0, 0)
    ppgplot.pgline(periods, power)
    ppgplot.pglab("Period (d)", "Amplitude", "Lomb-Scargle: " + arg.objectname)
    bestPeriod = periods[numpy.argmax(power)]
    lc = ppgplot.pgqci()
    ls = ppgplot.pgqls()
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([bestPeriod, bestPeriod], [0, max(power)])
    ppgplot.pgsci(lc)
    ppgplot.pgsls(ls)
    print "Best period: %f days or %f hours"%(bestPeriod, bestPeriod * 24.)
    

    


if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description='Loads a CSV file containing HJDs and RVs and tries to fit a period and sinusoid to the data.')
    # parser.add_argument('inputfile', type=str, help='Filename of the CSV file containing the RVs.')
    parser.add_argument('--device', type=str, default = "/xs", help='[Optional] PGPLOT device. Defaults to "/xs".')
    parser.add_argument('--title', type=str, help='Title for the plot. Otherwise title will be generated from data in the .CSV file.')
    parser.add_argument('--plo', type=float, default= 0.01, help='Period (days) to start the search. Default is 0.01 days.') 
    parser.add_argument('--phi', type=float, default= 10.0, help='Period (days) to stop the search. Default is 10 days.') 
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
        if good==0: print bcolors.WARNING + "%f\t%f\t%f"%(d['HJD'], d['RV'], d['RV error']) + bcolors.ENDC 
        else: 
            print "%f\t%f\t%f"%(d['HJD'], d['RV'], d['RV error'])
            dates.append(d['HJD'])
            velocities.append(d['RV'])
            velErrors.append(d['RV error'])
            
    print len(dates)
    
    
    phasePGPlotWindow = ppgplot.pgopen(arg.device)  
    ppgplot.pgask(False)
    pgPlotTransform = [0, 1, 0, 0, 0, 1]
    
    ppgplot.pgslct(phasePGPlotWindow)   
    ppgplot.pgsci(1)
    xStart = 2457000
    xStart = dates[0]
    reducedDates = [d - xStart for d in dates]
    startDate = min(reducedDates)
    endDate = max(reducedDates)
    maxVel = max(velocities)
    minVel = min(velocities)
    velRange = max([maxVel, abs(minVel)])
    print "Start date: %f, End date: %f, Velocity range: %f km/s"%(startDate, endDate, velRange)
    ppgplot.pgenv(0, len(dates), - 1.2 * velRange, 1.2 * velRange, 0, 0 )
    #ppgplot.pgenv(startDate, endDate, -velRange, velRange, 0, 0)
    
    ppgplot.pgpt(range(len(dates)), velocities)
    ppgplot.pgerrb(2, range(len(dates)), velocities, velErrors, 0)
    ppgplot.pgerrb(4, range(len(dates)), velocities, velErrors, 0)

    ppgplot.pglab("Point number", "Radial velocity km/s", arg.objectname)
    
    
    x = numpy.array(reducedDates)
    y = numpy.array(velocities)
    # Subtract the mean from the y-data
    y_mean = numpy.mean(y)
    y = y - y_mean
    import scipy.signal as signal
    
    periods = numpy.linspace(plo, phi, 1000)
    ang_freqs = 2 * numpy.pi / periods
    power = signal.lombscargle(x, y, ang_freqs)
    # normalize the power
    N = len(x)
    power *= 2 / (N * y.std() ** 2)

    pgramPGPlotWindow = ppgplot.pgopen(arg.device)  
    ppgplot.pgask(False)
    pgPlotTransform = [0, 1, 0, 0, 0, 1]

    ppgplot.pgslct(pgramPGPlotWindow)
    ppgplot.pgenv(min(periods), max(periods), 0, max(power), 0, 0)
    ppgplot.pgline(periods, power)
    ppgplot.pglab("Period (d)", "Amplitude", "Lomb-Scargle: " + arg.objectname)
    bestPeriod = periods[numpy.argmax(power)]
    lc = ppgplot.pgqci()
    ls = ppgplot.pgqls()
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([bestPeriod, bestPeriod], [0, max(power)])
    ppgplot.pgsci(lc)
    ppgplot.pgsls(ls)
    print "Best period: %f days or %f hours"%(bestPeriod, bestPeriod * 24.)
    
    # Now plot a folded RV curve
    t0 = reducedDates[0]
    periodDays = bestPeriod
    phaseFirst = [((d - t0) % periodDays) / periodDays for d in reducedDates]
    velocityPlot = copy.deepcopy(velocities)
    velocityErrorPlot = copy.deepcopy(velErrors)
    phases = copy.deepcopy(phaseFirst)
    for index, p in enumerate(phaseFirst):
        phases.append(p + 1.0)
        velocityPlot.append(velocities[index])
        velocityErrorPlot.append(velErrors[index])
    
    ppgplot.pgslct(phasePGPlotWindow)
    ppgplot.pgsci(1)
    ppgplot.pgenv(0, 2.0, -velRange*1.2, velRange*1.2, 0, 0)
    ppgplot.pgpt(phases, velocityPlot)
    ppgplot.pgerrb(2, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pgerrb(4, phases, velocityPlot, velocityErrorPlot, 0)
    ppgplot.pglab("Phase", "Radial velocity km/s", arg.objectname)
    
    # Now fit a sine wave with this period to the data...
    frequency = 1 / periodDays
    phase = 0
    gamma = 0
    amplitude = velRange
    
    x_values = numpy.array(phaseFirst)
    y_values = numpy.array(velocities)
    y_errors = numpy.array(velErrors)
    
    guess = [gamma, amplitude, phase]
    upperBounds = [100, 200, 1]
    lowerBounds = [-100, 0 , 0]
    bounds = (lowerBounds, upperBounds)
    print "Guess: ", guess
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
    
    loop = True
    while loop:
        print "Choice: [A] - Accept current period for least squares fit."
        print "\t[Z] - Zoom in the periodogram."
        print "\t[X] - Exit."
        choice = sys.stdin.read(1)
        if choice == 'x': sys.exit()
        if choice == 'z':
            print "Choose a period to zoom in by clicking on the periodogram."
            ppgplot.pgslct(pgramPGPlotWindow)
            (x1, y, char) = ppgplot.pgcurs(0, 0)
            print "Zooming in on range %f days to ...."%x1, 
            (x2, y, char) = ppgplot.pgband(4, 0, x1, y)
            range = phi-plo
            print x2, "days."
            getPeriodogram(pgramPGPlotWindow, reducedDates, velocities, x1, x2)
    
    ppgplot.pgslct(pgramPGPlotWindow)
    (x, y, char) = ppgplot.pgcurs(0, 0)
    print "Char: ", char
    if char=='q':
        periodGuess = bestPeriod
    else:
        periodGuess = x
    
    # Now try to tweak the period with a fixed amplitude and gamma velocity
    phaseGuess = 0
    guess = [periodGuess, phaseGuess]
    print "Period starting point is:", periodGuess
    gamma = gammaFit
    amplitude = amplitudeFit
    results, covariance = scipy.optimize.curve_fit(sineFixedGammaAmplitude, reducedDates, y_values, p0 = guess, sigma = y_errors)
    errors = numpy.sqrt(numpy.diag(covariance))
    
    tweakedPeriod = results[0]
    tweakedPeriodError = errors[0]
    phase = results[1]
    phaseError = errors[1]
    print "Tweaked Period: %f[%f]"%(tweakedPeriod, tweakedPeriodError)
    print "Tweaked Phase: %f[%f]"%(phase, phaseError)
    
    phases = [(d % tweakedPeriod)/tweakedPeriod for d in reducedDates]
    extendedPhases = numpy.asarray(phases)
    extendedVelocities = numpy.asarray(y_values)
    extendedVelocityErrors = numpy.asarray(y_errors)
    for index, p in enumerate(phases):
        extendedPhases = numpy.append(extendedPhases, p+1)
        extendedVelocities = numpy.append(extendedVelocities, y_values[index])
        extendedVelocityErrors = numpy.append(extendedVelocityErrors, y_errors[index])
    
    ppgplot.pgslct(mainPGPlotWindow)
    ppgplot.pgask(False)
    ppgplot.pgenv(0, 2.0, -velRange*1.2, velRange*1.2, 0, 0)
    yMin = gamma - 1.2*amplitude
    yMax = gamma + 1.2*amplitude
    ppgplot.pgenv(0, 2.0, yMin, yMax, 0, 0)
    periodHours = tweakedPeriod * 24.
    periodHoursError = tweakedPeriodError * 24.
    ppgplot.pglab("Phase", "Radial velocity km/s", arg.objectname + " period: %s[%s]d"%generalUtils.formatValueError(tweakedPeriod, tweakedPeriodError) + " or %s[%s] hours."%generalUtils.formatValueError(periodHours, periodHoursError))
    ppgplot.pgpt(extendedPhases, extendedVelocities, 0)
    ppgplot.pgerrb(2, extendedPhases, extendedVelocities, extendedVelocityErrors, 0)
    ppgplot.pgerrb(4, extendedPhases, extendedVelocities, extendedVelocityErrors, 0)
    ppgplot.pgsci(3)
    ppgplot.pgsls(2)
    ppgplot.pgline([0, 2], [gammaFit, gammaFit])
    ppgplot.pgsci(1)
    ppgplot.pgsls(4)
    ppgplot.pgline([0, 2], [0, 0])
  
    xFit = numpy.arange(0, 2, 0.01)
    yFit = gamma + amplitude * numpy.sin(numpy.pi * 2.0 * xFit + phase)
    x0 = - phase * 2 * numpy.pi / tweakedPeriod
    print "x0", x0
    print reducedDates
    print dates
    t0 = x0 + xStart
    print "T0", t0
    ppgplot.pgsci(2)
    ppgplot.pgsls(1)
    ppgplot.pgline(xFit, yFit)
    
    
    if not generalUtils.query_yes_no("Are you happy with this period?", default="no"):
        
        sys.exit()
    
 
    ppgplot.pgend()
    sys.exit()
    
   
