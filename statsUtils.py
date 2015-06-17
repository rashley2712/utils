import math
import numpy



def binDataWithErrors(x, y, yErrors, binFactor):
	newX = []
	newY = []
	newErrors = []
	
	newLength = len(x)/binFactor
	
	for i in range(newLength):
		xSum = 0
		ySum = 0
		weightSum = 0
		
		for j in range(binFactor):
			xSum+= x[i*binFactor + j]
			yVal = y[i*binFactor + j]
			yError = yErrors[i*binFactor + j]
			weight = 1/(yError*yError)
			ySum+= weight * yVal
			weightSum+= weight
			
		xVal = xSum/binFactor
		yVal = ySum/weightSum
		yError = math.sqrt(1/weightSum)
		newX.append(xVal)
		newY.append(yVal)
		newErrors.append(yError)
	
	return newX, newY, newErrors
	
def findNearestTime(time, timeValues, yValues):
	distance = 1000
	closestIndex = 0
	for index, t in enumerate(timeValues):
		gap = abs(time - t)
		if gap < distance:
			distance = gap
			closestIndex = index
		
	return closestIndex, distance

