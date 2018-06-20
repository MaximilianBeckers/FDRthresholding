import numpy as np
import subprocess
import math
import gc
import os

#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group

#-------------------------------------------------------------------------------------
def estimateNoiseFromMap(map, windowSize, boxCoord):

	#**************************************************
	#****** function to estimate var an mean from *****
	#**** nonoverlapping boxes outside the particle ***
	#**************************************************

	if boxCoord == 0:
		#extract a sample of pure noise from the map
		sizeMap = map.shape;
		sizePatch = np.array([windowSize, windowSize, windowSize]);
		center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]]);
		sampleMap1 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(0.02*sizeMap[1]):(int(0.02*sizeMap[1]) + sizePatch[1]),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];

		sampleMap2 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(0.98*sizeMap[1] - sizePatch[1]):(int(0.98*sizeMap[1])),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];
	
		sampleMap3 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
		int(0.02*sizeMap[2]):(int(0.02*sizeMap[2]) + sizePatch[2])];

		sampleMap4 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
		int(0.98*sizeMap[2]) - sizePatch[2]:(int(0.98*sizeMap[2]))];

		#concatenate the two samples
		sampleMap = np.concatenate((sampleMap1, sampleMap2, sampleMap3, sampleMap4), axis=0);

	else:
		sizePatch = np.array([windowSize, windowSize, windowSize]);
		center = np.array(boxCoord);
		sampleMap = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(center[1]-0.5*sizePatch[1]):(int(center[1]-0.5*sizePatch[1]) + sizePatch[1]),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];	

	
	#estimate variance and mean from the sample
	mean = np.mean(sampleMap);
	var = np.var(sampleMap);
	
	return mean, var, sampleMap;
#-------------------------------------------------------------------------------------
def estimateNoiseFromMapInsideMask(map, mask):

	#**************************************************
	#****** function to estimate var an mean from *****
	#******* map outside the user provided mask *******
	#**************************************************

	mask[mask<=0.0] = 0.0;
	mask[mask>0.0] = 1000.0;
	mask[mask<1000.0] = 1.0;
	mask[mask==1000.0] = 0.0;

	sampleMap = np.copy(map)*mask;
	sampleMap = sampleMap[sampleMap != 0.0];
	
	#estimate variance and mean from the sample
	mean = np.mean(sampleMap);
	var = np.var(sampleMap);
	
	
	return mean, var, sampleMap;
#-----------------------------------------------------------------------------------
def estimateNoiseFromMapTomo(map):

	#**************************************************
	#************ estimate noise from tomogram ********
	#**************************************************

	#estimate noise from whole tomogram by means of robust statistics
	median = np.median(map);
	medianAbsoluteDeviation = np.median(np.abs(map - median));
	
	#scale the median to get consistent estimator for the standard deviation
	sdRobust = 1.4826 * medianAbsoluteDeviation;
	varRobust = sdRobust * sdRobust;

	#varRobust = np.var(map);
	return median, varRobust;

#-----------------------------------------------------------------------------------
def estimateNoiseFromMap2D(map, windowSize):

	#**************************************************
	#****** function to estimate var an mean from *****
	#**** nonoverlapping boxes outside the particle ***
	#**************************************************

	#extract a sample of pure noise from the map
	sizeMap = map.shape;
	sizePatch = np.array([windowSize, windowSize]);
	center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1]]);
	sampleMap1 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
	int(0.02*sizeMap[1]):(int(0.02*sizeMap[1]) + sizePatch[1])];

	sampleMap2 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
	int(0.98*sizeMap[1] - sizePatch[1]):(int(0.98*sizeMap[1]))];

	#concatenate the two samples
	sampleMap = np.concatenate((sampleMap1, sampleMap2), axis=0);

	#estimate white noise in the map
	mean = np.mean(sampleMap);
	#print('Estimated mean of the background noise: {}'.format(mean))
	var = np.var(sampleMap);
	#print('Estimated variance of the mean background noise: {}'.format(var))

	return mean, var;
#-------------------------------------------------------------------------------------
def estimateECDFFromMap(map, windowSize, boxCoord):

	#**************************************************
	#****** function to estimate empirical cumul. *****
	#**** distribution function from solvent area *****
	#**************************************************

	if boxCoord == 0:
		#extract a sample of pure noise from the map
		sizeMap = map.shape;
		sizePatch = np.array([windowSize, windowSize, windowSize]);
		center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]]);
		sampleMap1 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(0.02*sizeMap[1]):(int(0.02*sizeMap[1]) + sizePatch[1]),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];

		sampleMap2 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(0.98*sizeMap[1] - sizePatch[1]):(int(0.98*sizeMap[1])),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];
	
		sampleMap3 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
		int(0.02*sizeMap[2]):(int(0.02*sizeMap[2]) + sizePatch[2])];

		sampleMap4 = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
		int(0.98*sizeMap[2]) - sizePatch[2]:(int(0.98*sizeMap[2]))];

		#conatenate the two samples
		sampleMap = np.concatenate((sampleMap1, sampleMap2, sampleMap3, sampleMap4), axis=0);
	
	elif boxCoord == -1:
		sampleMap = map;
	
	else:
		sizePatch = np.array([windowSize, windowSize, windowSize]);
		center = np.array(boxCoord);
		sampleMap = map[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
		int(center[1]-0.5*sizePatch[1]):(int(center[1]-0.5*sizePatch[1]) + sizePatch[1]),
		(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))];	
	
	#estimate ECDF from map
	sampleMap = sampleMap.flatten();	
	
	#downsize the sample
	finalSampleSize = min(100000, sampleMap.size);
	sampleMap = np.random.choice(sampleMap, finalSampleSize);
	numSamples = sampleMap.size;
	sampleMapSort = np.sort(sampleMap);

	minX = sampleMapSort[0];
	maxX = sampleMapSort[numSamples-1];
	numInterval = numSamples;
	spacingX = (maxX - minX)/(float(numInterval));

	ECDF = np.zeros(numInterval);	
	for index in range(numInterval):
		#val = minX + spacingX*index;
		val = sampleMapSort[index];
		ECDF[index] = ((sampleMapSort[sampleMapSort<= val]).size)/float(numSamples);
	
	return ECDF, sampleMapSort;

#------------------------------------------------------------------------------------
def getCDF(x, ECDF, sampleMapSort):

	#****************************************************
	#********* get the value of the CDF at point x ******
	#******* CDF : Cumulative distribution function *****
	#****************************************************
	
	numSamples = sampleMapSort.size;
	minX = sampleMapSort[0];
	maxX = sampleMapSort[numSamples-1];
	numInterval = 10000;

	spacingX = (maxX - minX)/(float(numInterval));

	if x >= maxX:
		CDFval = 1.0;
	elif x <= minX:
		CDFval = 0.0;
	else:
		#get the index in the ECDF array		
		index = int(math.floor((x - minX)/float(spacingX)));	
		CDFval = ECDF[index];

	return CDFval;

#------------------------------------------------------------------------------------
def normalizeMap(map, mean, var):
	
	#****************************************
	#********* normalize map ****************
	#****************************************
	
	if np.isscalar(var):
	 	normMap = np.subtract(map, mean);
		normMap = np.multiply(normMap, (1.0/(math.sqrt(var))));
	else: #if local variances are known, use them
		var[var==0] = 1000;
		normMap = np.subtract(map, mean);
		normMap = np.divide(normMap, np.sqrt(var));

	return normMap;


#-----------------------------------------------------------------------------------
def calcQMap(map, mean, var, ECDF, windowSize, boxCoord, mask, method, test):

	#*****************************************
	#***** generate qMap of a 3D density *****
	#*****************************************

	#get some map data
	sizeMap = map.shape;

	#calculate the test statistic
	if np.isscalar(var):
		map[map == 0.0] = -100000000;
		map = np.subtract(map, mean);   
		tMap = np.multiply(map, (1.0/(math.sqrt(var))));
		map = np.copy(tMap);
	else:
		var[var==0.0] = 1000.0; #just to avoid division by zero
		map = np.subtract(map, mean); 	
		tMap = np.divide(map, np.sqrt(var));
		var[var==1000.0] = 0.0;
		#upadte the mask, necessary as resmap is masking as well
		mask = np.multiply(mask,var);

	if test == 'rightSided':
		tMap[tMap<0] = -10000000.0;
	elif test == 'leftSided':
		tMap[tMap>0] = 10000000.0;

	#calculate the p-Values
	print('Calculating p-Values ...');
	pMap = np.zeros(sizeMap);
	
	if np.isscalar(ECDF): 
		if ECDF == 0: 
			vectorizedErf = np.vectorize(math.erf);
			erfMap = vectorizedErf(tMap/math.sqrt(2.0));
			#erf2Map = special.erf(tMap/math.sqrt(2.0));

			pMapRight = 1.0 - (0.5*(1.0 + erfMap)); 
			pMapLeft = (0.5*(1.0 + erfMap));
	
		else:
			#if ecdf shall be used, use if to p-vals
			ECDF, sampleSort = estimateECDFFromMap(map, windowSize, boxCoord);
			print('start ECDF calculation ...');	
			vecECDF = np.interp(map, sampleSort, ECDF, left=0.0, right=1.0);
			pMapRight = 1.0 - vecECDF;
			pMapLeft = vecECDF;

			#do test plots
			X = np.linspace(-5, 5, 10000);
			vectorizedErf = np.vectorize(math.erf);
			Y1 = np.interp(X, sampleSort, ECDF, left=0.0, right=1.0);
			Y2 = (0.5*(1.0 + vectorizedErf(X/math.sqrt(2.0))));
			import matplotlib.pyplot as plt
			plt.plot(X,Y1);
			plt.plot(X,Y2);
			plt.savefig('ECDF_vs_CDF.png', dpi=600);
			plt.close();
	    	#do test plots
			X = np.linspace(2, 5, 10000);
			Y1 = np.interp(X, sampleSort, ECDF, left=0.0, right=1.0);
			Y2 = (0.5*(1.0 + vectorizedErf(X/math.sqrt(2.0))));
			import matplotlib.pyplot as plt 
			plt.plot(X,Y1);
			plt.plot(X,Y2);
			plt.savefig('ECDF_vs_CDF_tail.png', dpi=600);
			plt.close();
	else:
		pMapRight = 1 - ECDF;
		pMapLeft = ECDF;


	if test == 'twoSided':
		pMap = np.minimum(pMapLeft, pMapRight);
	elif test == 'rightSided':
		pMap = pMapRight;
		pMap[tMap==-10000000.0] = 1.0
	elif test == 'leftSided':
		pMap = pMapLeft;
		pMap[tMap==10000000.0] = 1.0;
	
	#take the p-values in the mask	
	binaryMask = np.copy(mask);
	binaryMask[ binaryMask != 0.0 ] = 1.0;
	binaryMask[ binaryMask ==  0.0] = np.nan;

	pMap = pMap * binaryMask;
	pFlat = pMap.flatten();
	pInBall = pFlat[pFlat != np.nan];

	#do FDR control, i.e. calculate the qMap
	print('Start FDR control ...');
	qValues = FDR(pInBall, method);	

	qFlat = np.copy(pFlat);
	qFlat[qFlat != np.nan] = qValues;
	qMap = np.reshape(qFlat, sizeMap);
	qMap[qMap == np.nan] = 1.0;

	return qMap;

#---------------------------------------------------------------------------------
def FDR(pValues, method):

	#***********************************
	#***** FDR correction methods ******
	#***********************************

	numPVal = len(pValues);

	pSortInd = np.argsort(pValues);
	pSort = pValues[pSortInd];
	
	pAdjust = np.zeros(numPVal);
	prevPVal = 1.0;

	#use expansion for harmonic series
	Hn = math.log(numPVal) + 0.5772 + 0.5/numPVal - 1.0/(12*numPVal**2) + 1.0/(120*numPVal**4);    

	if method =='BH': #do benjamini-hochberg procedure
		for i in range(numPVal-1, -1, -1):
			pAdjust[i] = min(prevPVal, pSort[i]*numPVal/(i+1.0));
			prevPVal = pAdjust[i];	
	
	elif method == 'BY': #do benjamini yekutieli procedure
		for i in range(numPVal-1, -1, -1):
			pAdjust[i] =  min(prevPVal, pSort[i]*(numPVal/(i+1.0))*Hn);
			prevPVal = pAdjust[i];
	else:
		print('Please specify a method');
		return;	

	#sort back to the original order
	pSortIndOrig = np.argsort(pSortInd);

	return pAdjust[pSortIndOrig];

#---------------------------------------------------------------------------------
def binarizeMap(map, threshold):

	#***********************************
	#*binarize map at given threshold **
	#***********************************

	binMap = np.array(map);
	binMap[binMap <= threshold] = 0;
	binMap[binMap > threshold] = 1;

	#finally invert the map
	binMap = np.subtract(np.ones(binMap.shape), binMap);
  	
	return binMap;

#---------------------------------------------------------------------------------
def printSummary(args, time):

	#***********************************
	#** print a Summary of the job *****
	#***********************************

	print("*** Done ***");
	print(" ");
	print("******** Summary ********");

	output = "Elapsed Time: " + repr(time);
	print(output);

	#print input map filename
	splitFilename = os.path.splitext(os.path.basename(args.em_map));
	output = "Input EM-map: " + splitFilename[0] + ".mrc"; 
	print(output);

	#print model map filename
	if args.model_map is not None:
		splitFilename = os.path.splitext(os.path.basename(args.model_map));
		output = "LocScale was done with the input model-map: " + splitFilename[0] + ".mrc";
		print(output);

		#print window size
		if args.window_size is not None:
			w = int(math.ceil(args.window_size / 2.) * 2);
		else:
			w = int(round(7 * 3 * args.apix));
		output = "Window size: " + repr(w);

	#print local resolution map filename
	if args.locResMap is not None:
		splitFilename = os.path.splitext(os.path.basename(args.locResMap));
		output = "Input local resolution map: " + splitFilename[0] + ".mrc";
		print(output);
	
	#print output filenames
	if args.outputFilename is not None:
		splitFilename = os.path.splitext(os.path.basename(args.outputFilename));
	else:
		splitFilename = os.path.splitext(os.path.basename(args.em_map));
	
	
	if args.model_map is not None:
		output = "Output LocScale map: " + splitFilename[0] + "_scaled" + ".mrc";
		print(output);
		output = "Output confidence Map: " + splitFilename[0] + "_confidenceMap" + ".mrc";
		print(output);
	else:
		if args.locResMap is not None:
			output = "Output locally filtered map: " + splitFilename[0] + "_locFilt" + ".mrc";
			print(output);
			output = "Output confidence Map: " + splitFilename[0] + "_confidenceMap" + ".mrc";
			print(output);
		else:
			output = "Output confidence Map: " + splitFilename[0] + "_confidenceMap" + ".mrc";
			print(output);

	#print pixel size
	output = "Pixel size: " + repr(args.apix);
	print(output);
	
				


