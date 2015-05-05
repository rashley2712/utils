import matplotlib.pyplot
import ppgplot, numpy 

def matplot(slot):
	print "We are about to plot: %s"%(slot)
	if plotterHandle!=None:
		figure = plotterHandle
		matplotlib.pyplot.figure(figure)
	else:
		figure = matplotlib.pyplot.figure(figsize=(12, 8))
	
	xValues = slot.photometry.times
	yValues = slot.getColumn("Counts")
	
	matplotlib.pyplot.plot(xValues, yValues, 'r.')
	matplotlib.pyplot.ioff()
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block=False)
	return figure
	
def pgplot(slot):
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
	