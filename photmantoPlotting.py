import matplotlib.pyplot
import ppgplot, numpy 

def matplot(slot, state):
	print "We are about to plot, %s using matplotlib"%(slot)
	
	if state['plotterhandle']!=None:
		figNumber = state['plotterhandle']
		matplotlib.pyplot.figure(figNumber)
	else:
		print "Creating a new figure."
		figure = matplotlib.pyplot.figure(figsize=(12, 8))
		state['plotterhandle'] = figure.number
	
	if not state['overplot']:
		matplotlib.pyplot.clf() 		# Clear the current plot
	
	xValues = slot.photometry.times
	yValues = slot.getColumn("Counts")
	
	plotSymbols = state['plotcolour'] + '.'
	matplotlib.pyplot.plot(xValues, yValues, plotSymbols)
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block=False)
	return state
	
def pgplot(slot, plotterHandle):
	print "We are about to plot: %s"%(slot)
	
	lightcurveView = ppgplot.pgopen('/xs')
	
	xValues = slot.photometry.times
	xValues = numpy.arange(0, len(xValues))
	yValues = slot.getColumn("Counts")
	
	ppgplot.pgenv(min(xValues), max(xValues), min(yValues), max(yValues), 0, 0)
	ppgplot.pgask(False)
	
	ppgplot.pgsci(2)	
	ppgplot.pgpt(xValues, yValues, 1)

	ppgplot.pgclos()
	return
	