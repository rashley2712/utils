import matplotlib.pyplot
import ppgplot

def matplotlibplot(slot):
	print "We are about to plot: %s"%(slot)
	
	figure = matplotlib.pyplot.figure(figsize=(12, 8))
	
	xValues = slot.photometry.times
	yValues = slot.getColumn("Counts")
	
	matplotlib.pyplot.plot(xValues, yValues, 'r.')
	
	matplotlib.pyplot.show()
	return
	
def pgplot(slot):
	print "We are about to plot: %s"%(slot)
	
	lightcurveView = ppgplot.pgopen('/xs')
	
	xValues = slot.photometry.times
	yValues = slot.getColumn("Counts")
	
	ppgplot.pgenv(min(xValues), max(xValues), min(yValues), max(yValues), 0, 0)
	ppgplot.pgask(False)
	
	ppgplot.pgsci(2)	
	ppgplot.pgpt(xValues, yValues, 1)

	ppgplot.pgclos()
	return
	