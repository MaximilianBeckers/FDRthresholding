from EMAN2 import *
from FDRutil import *
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import argparse, os, sys
import subprocess
import math
import gc
import os.path
from time import sleep


#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group (2017)


#------------------------------------------------------------------------------------------------------
def localFiltration(map, locResMap, apix, localVariance, windowSize):

	#**************************************************
	#**** function to perform a local filtration ******
	#****** according to the local resolution *********
	#**************************************************


	#some initialization
	numX, numY, numZ = map.get_xsize(), map.get_ysize(), map.get_zsize();
	mean = np.zeros((numX, numY, numZ));
	var = np.zeros((numX, numY, numZ));
	filteredMapData = np.zeros((numX, numY, numZ));

	#transform to numpy array
	locResMapData = EMNumPy.em2numpy(locResMap);

	#transform to abosulte frequency units(see http://sparx-em.org/sparxwiki/absolute_frequency_units)
	locResMapData = np.divide(apix, locResMapData);
	
	#round to 3 decimals
	locResMapData = np.around(locResMapData, 3);	

	#set resolution search range, 3 decimals exact
	locResArray = np.arange(0, 0.5 , 0.001);

	counter = 0;
	numRes = len(locResArray);	

	#get initial noise statistics
	initMapData = np.copy(EMNumPy.em2numpy(map));
	initMean, initVar, _ = estimateNoiseFromMap(initMapData, windowSize);
	noiseMapData = np.random.normal(initMean, math.sqrt(initVar), (100, 100, 100));
	noiseMap = EMNumPy.numpy2em(noiseMapData);

	# Initial call to print 0% progress
	#printProgressBar(counter, numRes, prefix = 'Progress:', suffix = 'Complete', bar_length = 50)
	
	for tmpRes in locResArray:   
		
		#get indices of voxels with the current resolution	
		indices = np.where(locResMapData == tmpRes);
	
		if (indices[0].size == 0):
			#this resolution is obviously not in the map, so skip
			counter = counter + 1;
			continue;
		elif (tmpRes == round(apix/100.0, 3)):
			xInd, yInd, zInd = indices[0], indices[1], indices[2];
			
			#do local filtration
			tmpFilteredMap = map.process("filter.lowpass.tanh", {"cutoff_abs": tmpRes, "fall_off": 0.1});
			tmpFilteredMapData = np.copy(EMNumPy.em2numpy(tmpFilteredMap));

			#set the filtered voxels
			filteredMapData[xInd, yInd, zInd] = tmpFilteredMapData[xInd, yInd, zInd];

		else:
		
			xInd, yInd, zInd = indices[0], indices[1], indices[2];
			
			#do local filtration
			tmpFilteredMap = map.process("filter.lowpass.tanh", {"cutoff_abs": tmpRes, "fall_off": 0.1});
			tmpFilteredMapData = np.copy(EMNumPy.em2numpy(tmpFilteredMap));
			
			#set the filtered voxels
			filteredMapData[xInd, yInd, zInd] = tmpFilteredMapData[xInd, yInd, zInd];
			
			if localVariance == True:
				#estimate and set noise statistic
				tmpMean, tmpVar, _ = estimateNoiseFromMap(tmpFilteredMapData, windowSize);
			
				mean[xInd, yInd, zInd] = tmpMean;
				var[xInd, yInd, zInd] = tmpVar;

	return filteredMapData, mean, var;

#---------------------------------------------------------------------------------------------------
def padFourier(map, apix, finalPix):
	
	sizeMap = np.array([map.get_xsize(), map.get_ysize(), map.get_zsize()]);
	
	#size of the map
	lengthMap = sizeMap[1] * apix;
	finalSize = int(lengthMap/float(finalPix));
	diffSize = finalSize - sizeMap[1];

	if diffSize % 2 != 0 :
		finalSize = finalSize + 1;
		diffSize = finalSize - sizeMap[1];
	
	newPixSize = lengthMap/float(finalSize);
	print('The new pixelSize is: {}'.format(newPixSize));
	
	halfDiffSize = int(diffSize/2.0);

	#do fourier padding
	mapData = np.copy(EMNumPy.em2numpy(map));
	fftMap = np.fft.fftn(mapData);
	fftMapshift = np.fft.fftshift(fftMap);
	
	fftMapPad = np.zeros((finalSize, finalSize, finalSize), dtype=np.complex_);
	fftMapPad[halfDiffSize:(halfDiffSize + sizeMap[1]), halfDiffSize:(halfDiffSize + sizeMap[1]), halfDiffSize:(halfDiffSize + sizeMap[1])] = np.copy(fftMapshift);
	mapPadData = np.fft.ifftn(np.fft.ifftshift(fftMapPad)).real;
	
	mapPad = EMNumPy.numpy2em(mapPadData);

	return mapPad;

#--------------------------------------------------------------------------------------------------
def makeDiagnosticPlot(map, windowSize, padded, singleBox, boxCoord):

	#*************************************************************
	#*** function to make diagnostic plot of noise estimation ****
	#*************************************************************

	print("Generating diagnostic plot of noise estimation. Please have a look in 'diag_image.pdf' that the molecule does not fall into the region used for background noise estimation.")

	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.gridspec as gridspec

	tmpMap = EMData();
	tmpMap = map;
	mapData = EMNumPy.em2numpy(tmpMap);

	sizeMap = mapData.shape;
	sizePatch = np.array([windowSize, windowSize, windowSize]);
	center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]]);
	visMap = np.copy(mapData);
	noiseLabel = (np.mean(visMap) + 5* np.var(visMap));
	

	#if coordinates are provided, do singleBox estimation
	if boxCoord != 0:
		singleBox = True;
	
	if singleBox == False:
		visMap[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
			int(0.02*sizeMap[1]+padded):(int(0.02*sizeMap[1]+padded) + sizePatch[1]),
			(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))] = noiseLabel;

		visMap[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
			int(0.98*sizeMap[1] - padded - sizePatch[1]):(int(0.98*sizeMap[1] - padded)),
			(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))] = noiseLabel;

		visMap[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
			(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
			int(0.02*sizeMap[2] + padded):(int(0.02*sizeMap[2] + padded) + sizePatch[2])] = noiseLabel;

		visMap[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
			(int(center[1]-0.5*sizePatch[1])):(int((center[1]-0.5*sizePatch[1]) + sizePatch[1])), 
			int(0.98*sizeMap[2] - padded - sizePatch[2]):(int(0.98*sizeMap[2]) - padded)] = noiseLabel;

	else:
		if boxCoord != 0:
			visMap[int(boxCoord[0]+padded-0.5*sizePatch[0]):(int(boxCoord[0]+padded-0.5*sizePatch[0]) + sizePatch[0]),
				int(boxCoord[1]+padded-0.5*sizePatch[1]):int((boxCoord[1]+padded-0.5*sizePatch[1]) + sizePatch[1]),
				(int(boxCoord[2]+padded-0.5*sizePatch[2])):(int((boxCoord[2]+padded-0.5*sizePatch[2]) + sizePatch[2]))] = noiseLabel;
		else:	
			visMap[int(center[0]-0.5*sizePatch[0]):(int(center[0]-0.5*sizePatch[0]) + sizePatch[0]),
				int(0.02*sizeMap[1]+padded):(int(0.02*sizeMap[1]+padded) + sizePatch[1]),
				(int(center[2]-0.5*sizePatch[2])):(int((center[2]-0.5*sizePatch[2]) + sizePatch[2]))] = noiseLabel;
	
	if boxCoord == 0:
		sliceMapYZ = visMap[int(sizeMap[0]/2.0), :, :];
		sliceMapXZ = visMap[:, int(sizeMap[1]/2.0), :];
		sliceMapXY = visMap[:, :, int(sizeMap[2]/2.0)];
	else:
		sliceMapYZ = visMap[boxCoord[0], :, :];
		sliceMapXZ = visMap[:, boxCoord[1], :];
		sliceMapXY = visMap[:, :, boxCoord[2]];

	#make diagnostics plot	
	plt.gray(); #make grayscale images
	plt.rc('xtick', labelsize=8);    # fontsize of the tick labels
	plt.rc('ytick', labelsize=8);    # fontsize of the tick labels
	pp = PdfPages('diag_image.pdf');
	gs = gridspec.GridSpec(1, 3);
	
	ax1 = plt.subplot(gs[2]);
	ax1.set_title('Y-Z slice');
	ax1.imshow(sliceMapYZ);
	
	ax2 = plt.subplot(gs[1]);
	ax2.set_title('X-Z slice');
	ax2.imshow(sliceMapXZ);

	ax3 = plt.subplot(gs[0]);
	ax3.set_title('X-Y slice');
	ax3.imshow(sliceMapXY);

	pp.savefig();
	pp.close();


#------------------------------------------------------------------------------------------------
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, bar_length = 100):
     
    #******************************************
    #** progress bar for local visualization **
    #******************************************	
    """
    Call in a loop to create terminal progress bar
    params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    str_format = "{0:." + str(decimals) + "f}";
    percents = str_format.format(100 * (iteration / float(total)));
    filled_length = int(round(bar_length * iteration / float(total)));
    bar = '#' * filled_length + '-' * (bar_length - filled_length);

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix));

    if iteration == total:
        sys.stdout.write('\n');
    sys.stdout.flush();

#---------------------------------------------------------------------------------------------------
def sharpenMap(map, bFactor, resolutionCutoff):

	#******************************************
	#******* do usual bFactor sharpening ******
	#******************************************

	cutoffFreq = 1.0/resolutionCutoff;
	sharpenedMap = map.process("filter.lowpass.autob",{"cutoff_freq": cutoffFreq, "bfactor": bFactor}); 

	return sharpenedMap;

#--------------------------------------------------------------------------------------------------
def readAndFlattenImageStack(filename):
	
	#***********************************************
	#* read imageStack and return Nx(nX*nY) array **
	#***********************************************

	imageStack = EMData();
	imageStack.read_image(filename);
	nx, ny, numImages = imageStack.get_xsize(), imageStack.get_ysize(), imageStack.get_zsize();
	imageStackData = EMNumPy.em2numpy(imageStack);
	imageStack = []; #free memory

	#flatten each image and append to array of NxD, with N the number of images and D the size of the flat image
	imageStackData = np.reshape(imageStackData, (nx*ny, numImages), order='F' );		
	imageStackData = np.transpose(imageStackData);
	print(imageStackData.shape);

	return imageStackData;

#---------------------------------------------------------------------------------
def addNoiseToMap(map, varNoise):

	#*********************************
	#****** add Noise To Map *********
	#*********************************

	mapData = np.copy(EMNumPy.em2numpy(map));
	noiseMap = np.random.randn(map.get_xsize(), map.get_ysize(), map.get_xsize())*math.sqrt(varNoise);
	mapData = mapData + noiseMap;
	map = EMNumPy.numpy2em(np.copy(mapData));

	return map 

