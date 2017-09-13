import numpy as np
import subprocess
import math
import gc
import os

#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group

#-------------------------------------------------------------------------------------
def estimateNoiseFromMap(map, windowSize):

	#**************************************************
	#****** function to estimate var an mean from *****
	#**** nonoverlapping boxes outside the particle ***
	#**************************************************

	#extract a sample of pure noise from the map
	sizeMap = map.shape
	sizePatch = np.array([windowSize, windowSize, windowSize])
	center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]])
	sampleMap1 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
	int(0.02*sizeMap[1]):(int(0.02*sizeMap[1]) + sizePatch[1]),
	(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))]

	sampleMap2 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
	int(0.98*sizeMap[1] - sizePatch[1]):(int(0.98*sizeMap[1])),
	(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))]
	
	sampleMap3 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
	(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
	int(0.02*sizeMap[2]):(int(0.02*sizeMap[2]) + sizePatch[2])]

	sampleMap4 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
	(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
	int(0.98*sizeMap[2]) - sizePatch[2]:(int(0.98*sizeMap[2]))]

	#concatenate the two samples
	sampleMap = np.concatenate((sampleMap1, sampleMap2, sampleMap3, sampleMap4), axis=0)

	#estimate white noise in the map
	mean = np.mean(sampleMap)
	#print('Estimated mean of the background noise: {}'.format(mean))
	var = np.var(sampleMap)
	#print('Estimated variance of the mean background noise: {}'.format(var))

	return mean, var

#-----------------------------------------------------------------------------------
def calcQMap(map, mean, var, mask, method, test):

	#*****************************************
	#***** generate qMap of a 3D density *****
	#*****************************************


	#get some map data
	sizeMap = map.shape

	#calculate the test statistic
	if np.isscalar(var):
		map = np.subtract(map, mean)    
		tMap = np.multiply(map, (1.0/(math.sqrt(var))))
	else:
		var[var==0] = 1000
		map = np.subtract(map, mean) 	
		tMap = np.divide(map, np.sqrt(var))

	#calculate the p-Values
	print('Calculating p-Values ...')
	pMap = np.zeros(sizeMap)

	pFlat = pMap.flatten()

	pValueCounter = 0;
	for iInd in range(0,sizeMap[0]):
		for jInd in range(0,sizeMap[1]):
			for kInd in range(0,sizeMap[2]):
				
				if (mask[iInd, jInd, kInd] != 0) :
					pValueRight = 1.0 - (0.5 *(1.0 +  math.erf(tMap[iInd, jInd, kInd]/math.sqrt(2.0))))	
					pValueLeft = (0.5 *(1.0 +  math.erf(tMap[iInd, jInd, kInd]/math.sqrt(2.0))))					

					if test == 'twoSided':
						pValue = min(pValue, pValueDown)

					elif test == 'rightSided':
						pValue = pValueRight
			
					elif test == 'leftSided':
						pValue = pValueLeft
		
					pFlat[pValueCounter] = pValue
				else:
					pFlat[pValueCounter] = -1
				
				pValueCounter = pValueCounter + 1

	pInBall = pFlat[pFlat != -1]
	
	if np.isscalar(var) == False :
		var[var==1000] = 0    

	#do FDR control, i.e. calculate the qMap
	print('Start FDR control ...')
	qValues = FDR(pInBall, method);	
                
	qMap = np.ones(sizeMap)			
	qValueCounter = 0
	pValueCounter = 0
	for iInd in range(0,sizeMap[0]):
		for jInd in range(0,sizeMap[1]):
			for kInd in range(0,sizeMap[2]):
				if pFlat[pValueCounter] != -1:
					qMap[iInd, jInd, kInd] = qValues[qValueCounter]
	 				qValueCounter = qValueCounter + 1               
				
				pValueCounter = pValueCounter + 1

	return qMap


#---------------------------------------------------------------------------------
def FDR(pValues, method):

	#***********************************
	#***** FDR correction methods ******
	#***********************************


	numPVal = len(pValues)

	pSortInd = np.argsort(pValues)
	pSort = pValues[pSortInd]
	
	pAdjust = np.zeros(numPVal)
	prevPVal = 1.0

	#use expansion for harmonic series
	Hn = math.log(numPVal) + 0.5772 + 0.5/numPVal - 1.0/(12*numPVal**2) + 1.0/(120*numPVal**4)    

	if method =='BH': #do benjamini-hochberg procedure
		for i in range(numPVal-1, -1, -1):
			pAdjust[i] = min(prevPVal, pSort[i]*numPVal/(i+1.0))
			prevPVal = pAdjust[i]	
	
	elif method == 'BY': #do benjamini yekutieli procedure
		for i in range(numPVal-1, -1, -1):
			pAdjust[i] =  min(prevPVal, pSort[i]*(numPVal/(i+1.0))*Hn)
			prevPVal = pAdjust[i]
	else:
		print('Please specify a method')
		return	

	#sort back to the original order
	pSortIndOrig = np.argsort(pSortInd)

	return pAdjust[pSortIndOrig]


#---------------------------------------------------------------------------------
def binarizeMap(map, threshold):

	#***********************************
	#*binarize map at given threshold **
	#***********************************

	binMap = np.array(map)
	binMap[binMap <= threshold] = 0
	binMap[binMap > threshold] = 1

	#finally invert the map
	binMap = np.subtract(np.ones(binMap.shape), binMap)
  	
	return binMap
	
#---------------------------------------------------------------------------------
def printSummary(args, time):

	#***********************************
	#** print a Summary of the job *****
	#***********************************

	print("*** Done ***")
	print(" ")
	print("******** Summary ********")


	output = "Elapsed Time: " + repr(time)
	print(output)

	#print input map filename
	splitFilename = os.path.splitext(os.path.basename(args.em_map))
	output = "Input EM-map: " + splitFilename[0] + ".mrc" 
	print(output)

	#print model map filename
	if args.model_map is not None:
		splitFilename = os.path.splitext(os.path.basename(args.model_map))
		output = "LocScale was done with the input model-map: " + splitFilename[0] + ".mrc"
		print(output)

		#print window size
		if args.window_size is not None:
			w = int(math.ceil(args.window_size / 2.) * 2)
		else:
			w = int(round(7 * 3 * args.apix))
		output = "Window size: " + repr(w)

	#print local resolution map filename
	if args.locResMap is not None:
		splitFilename = os.path.splitext(os.path.basename(args.locResMap))
		output = "Input local resolution map: " + splitFilename[0] + ".mrc"
		print(output)
	
	#print output filenames
	if args.model_map is not None:
		splitFilename = os.path.splitext(os.path.basename(args.em_map))
		output = "Output LocScale map: " + splitFilename[0] + "_locscale" + ".mrc"
		print(output)
		output = "Output LocScale qMap: " + splitFilename[0] + "_locscale_qMap_inv" + ".mrc"
		print(output)
	else:
		if args.locResMap is not None:
			splitFilename = os.path.splitext(os.path.basename(args.em_map))
			output = "Output locally filtered map: " + splitFilename[0] + "_localFilt" + ".mrc"
			print(output)
			output = "Output qMap: " + splitFilename[0] + "_qMap_inv" + ".mrc"
			print(output)
		else:
			splitFilename = os.path.splitext(os.path.basename(args.em_map))
			output = "Output qMap: " + splitFilename[0] + "_qMap_inv" + ".mrc"
			print(output)

	#print pixel size
	output = "Pixel size: " + repr(args.apix)
	print(output)
	
				


