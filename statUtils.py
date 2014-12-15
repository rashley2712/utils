import math
import numpy



def binDataWithErrors(x, y, yErrors, binFactor):
	newX = []
	newY = []
	newErrors = []
	
	newLength = len(x)/binFactor
	
	for i in range(newLength):
		xVal = 0
		yVal = 0
		for j in range(bin):
			xVal+= x[i*bin +j]
			yVal+= y[i*bin +j]
		xVal = xVal/bin
		yVal = yVal/bin
		newX.append(xVal)
		newY.append(yVal)
	
	return newX, newY
