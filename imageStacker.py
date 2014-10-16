#!/usr/bin/env python
import numpy, math
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse
import sys
from scipy.fftpack import fft2, ifft2
from scipy.ndimage.fourier import fourier_shift

def luminance(r, g, b):
	lum = math.sqrt( 0.299*r*r + 0.587*g*g + 0.114*b*b)
	return lum

def getLuminance(data):
	
	(x, y, depth) = numpy.shape(data)
	lumData = numpy.zeros((x, y))
	
	for i in range(x):
		for j in range(y):
			d = data[i][j]
			l = luminance(d[0], d[1], d[2])
			lumData[i][j] = l
		print i
		
	return lumData
	
def phase_cor( A, B ):
  """
  find the 2D correlation between images A and B. This is much faster than any
  of the other scipy correlation methods which insist on working in the spatial
  domain.
  """
  return ( ifft2( fft2(A)*numpy.conj(fft2(B)) ) ).real

def apply_shift( X, shift ):
  """
  Shift an image in the fourier domain
  """
  return ifft2( fourier_shift( fft2(X), shift ) ).real


def percentiles(data, lo, hi):
	""" Returns a normalised array where lo percent of the pixels are 0 and hi percent of the pixels are 255
	"""
	max = data.max()
	print "Max pixel value:", max
	dataArray = data.flatten()
	pHi = numpy.percentile(dataArray, hi)
	pLo = numpy.percentile(dataArray, lo)
	range = pHi - pLo
	scale = range/255
	data = numpy.clip(data, pLo, pHi)
	data-= pLo
	data/=scale
	return data

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Tries to stack images of starfields.')
	parser.add_argument('datafiles', nargs='+', type=str, help='Input data file(s)')
	
	arg = parser.parse_args()
	
	if len(arg.datafiles)<2 :
		print "I need at least two input files"
		sys.exit()
	
	colours = ['r', 'g', 'b']
	originalData = {'r': None, 'g': None, 'b': None}
	addedData = {'r': None, 'g': None, 'b': None}
	
	datafile = arg.datafiles[0]
	print "Opening the first file:", datafile
	imageData = matplotlib.pyplot.imread(datafile)
	print "Dimensions of the image:", numpy.shape(imageData)
	for i, c in enumerate(colours):
		print "Extracting colour:", c
		originalData[c] = imageData[:, :, i]
		addedData[c] = originalData[c]
				
	for datafile in arg.datafiles[1:]:
		print "Opening the next file:", datafile
		imageData = matplotlib.pyplot.imread(datafile)
		print "Dimensions of the image:", numpy.shape(imageData)
		for i, c in enumerate(colours):
			print "Extracting colour:", c, i
			newData = imageData[:, :, i]
			#shift = phase_cor(originalData[c], newData)
			# Find coordinates of the peak in the X-correlation
			#delta_calculated = numpy.unravel_index( numpy.argmax( shift ), numpy.shape( shift ) )
			#print "Shift measured to be", delta_calculated
			#shiftedData = apply_shift( newData, delta_calculated )
			shiftedData = newData
			addedData[c] = numpy.add(addedData[c],shiftedData)
			print addedData[c]
			
		for c in colours:
			print "Max value in %s is %d"%(c, numpy.max(addedData[c]))
		
	
	for c in colours:
		addedData[c] = float(addedData[c])/numpy.max(addedData[c]) * 255.
	
	for c in colours:
		print "Max value in %s is %d"%(c, numpy.max(addedData[c]))
		
	
	matplotlib.pyplot.figure()
	matplotlib.pyplot.title("Red image")
	image = matplotlib.pyplot.imshow(addedData['r'], cmap='Reds', interpolation='none')
	matplotlib.pyplot.show(block=False)
	
	finalImageData = numpy.dstack( (addedData['r'], addedData['g'], addedData['b']) ) 
	
	matplotlib.pyplot.figure()
	matplotlib.pyplot.title("Final image")
	image = matplotlib.pyplot.imshow(finalImageData)
	matplotlib.pyplot.show(block=True)
	
	
	
	
