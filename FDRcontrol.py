#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group (2017)

#import some stuff
from FDRutil import *
from mapUtil import *
from locscaleUtil import *
import numpy as np
import argparse, os, sys
import subprocess
import math
import gc
import os.path
import time
import sys
import mrcfile

#*************************************************************
#****************** Commandline input ************************
#*************************************************************

cmdl_parser = argparse.ArgumentParser(
               prog=sys.argv[0],
               description='*** Analyse density ***',
               formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=30), add_help=True);

cmdl_parser.add_argument('-em', '--em_map', metavar="em_map.mrc",  type=str, required=True,
                         help='Input filename EM map');
cmdl_parser.add_argument('-halfmap2', '--halfmap2', metavar="halfmap2.mrc",  type=str, required=False,
                         help='Input filename halfmap 2');
cmdl_parser.add_argument('-mask', '--mask',metavar="mask",  type=str,
                         help='Input filename mask');
cmdl_parser.add_argument('-p', '--apix', metavar="apix",  type=float, required=True,
                         help='pixel Size of input map');
cmdl_parser.add_argument('-fdr', '--fdr', metavar="fdr",  type=float, required=False,
                         help='False Discovery Rate');
cmdl_parser.add_argument('-locResMap', metavar="locResMap.mrc", type=str, required=False,
						help='Input local Resolution Map');
cmdl_parser.add_argument('-method', metavar="method", type=str, required=False,
						help="Method for multiple testing correction. 'BY' for Benjamini-Yekutieli, 'BH' for Benjamini-Hochberg or 'Holm' for Holm FWER control");
cmdl_parser.add_argument('-mm', '--model_map', metavar="model_map.mrc", type=str, required=False,
						help="Input model map for model based amplitude scaling");
cmdl_parser.add_argument('-w', '--window_size', metavar="windowSize", type=float, required=False,
						help="Input window size for local Amplitude scaling and background noise estimation");
cmdl_parser.add_argument('-mpi', action='store_true', default=False,
						help="Set this flag if MPI should be used for the local amplitude scaling");
cmdl_parser.add_argument('-o', '--outputFilename', metavar="output.mrc", type=str, required=False,
						help="Name of the output");
cmdl_parser.add_argument('-noiseBox', metavar="[x, y, z]", nargs='+', type=int, required=False,
							help="Box coordinates for noise estimation");
cmdl_parser.add_argument('-meanMap', '--meanMap', type=str, required=False,
                            help="3D map of noise means to be used for FDR control");
cmdl_parser.add_argument('-varianceMap', '--varianceMap', type=str, required=False,
                            help="3D map of noise variances to be used for FDR control");
cmdl_parser.add_argument('-testProc', '--testProc', type=str, required=False,
                            help="choose between right, left and two-sided testing");
cmdl_parser.add_argument('-lowPassFilter', '--lowPassFilter', type=float, required=False,
							help="Low-pass filter the map at the given resoultion prior to FDR control");
cmdl_parser.add_argument('-ecdf', action='store_true', default=False,
                        help="Set this flag if the empricical cumulative distribution function should be used instead of the standard normal distribution");
cmdl_parser.add_argument('-w_locscale', '--window_size_locscale', metavar="windowSize_locScale", type=float, required=False, 
						help="Input window size for local Amplitude scaling");
cmdl_parser.add_argument('-stepSize', '--stepSize', metavar="stepSize_locScale", type=int, required=False, 
						help="Voxels to skip for local amplitude scaling");

#************************************************************
#********************** main function ***********************
#************************************************************

def main():

	start = time.time();

	#get command line input
	args = cmdl_parser.parse_args();

	#if a model map is given, local amplitude scaling will be done
	if args.model_map is not None:
		
		if args.stepSize > args.window_size_locscale:
			print("Step Size cannot be bigger than the window_size. Job is killed ...")
			return;

		launch_amplitude_scaling(args);
	else:	
		#no ampltidue scaling will be done
		print('************************************************');
		print('******* Significance analysis of EM-Maps *******');
		print('************************************************');

		#load the maps
		if args.halfmap2 is not None:
			if args.em_map is None:
				print("One half map missing! Exit ...")
				sys.exit();
			else:
				#load the maps
				filename = args.em_map;
				map1 = mrcfile.open(args.em_map, mode='r+');
				halfMapData1 = np.copy(map1.data);

				map2 = mrcfile.open(args.halfmap2, mode='r+');
				halfMapData2 = np.copy(map2.data);

				mapData = (halfMapData1 + halfMapData2)*0.5;
				halfMapData1 = 0;
				halfMapData2 = 0;
                                
		else:
			#load single map
			filename = args.em_map;
			map = mrcfile.open(filename, mode='r+');
			mapData = np.copy(map.data);

		apix = args.apix;

		#get boxCoordinates
		if args.noiseBox is None:
			boxCoord = 0;
		else:
			boxCoord = args.noiseBox;
		
		#set test procdure
		if args.testProc is not None:
			testProc = args.testProc;
		else:
			testProc = 'rightSided';

		#set ECDF
		if args.ecdf:
			ECDF = 1;
		else:
			ECDF = 0;

		#set output filename
		if args.outputFilename is not None:
			splitFilename = os.path.splitext(os.path.basename(args.outputFilename));
		else:
			splitFilename = os.path.splitext(os.path.basename(filename));

		sizeMap = mapData.shape;

		if args.lowPassFilter is not None:
			frequencyMap = calculate_frequency_map(mapData);
			providedRes = apix/float(args.lowPassFilter);
			mapData = lowPassFilter(np.fft.rfftn(mapData), frequencyMap, providedRes, mapData.shape);

		#handle FDR correction procedure
		if args.method is not None:
			method = args.method;
		else:
			#default is Benjamini-Yekutieli
			method = 'BY';

		if args.window_size is not None:
			wn = args.window_size;
			wn = int(wn);
			if wn < 20:
				print("Provided window size is quite small. Please think about potential inaccuracies of your noise estimates!");
		else: 
			wn = max(int(0.05*sizeMap[0]),10);

		#generate a circular Mask
		sphere_radius = (np.min(sizeMap) // 2);
		circularMaskData = makeCircularMask( np.copy(mapData), sphere_radius);

		#if mask is provided, take it
		if args.mask is not None:	
			mask = mrcfile.open(args.mask, mode='r+');
			maskData = np.copy(mask.data);
		else:
			maskData = circularMaskData;
			
		#do some diagnostics and check for normality of map	
		makeDiagnosticPlot(mapData, wn, 0, False, boxCoord);
		checkNormality(mapData, wn, boxCoord);

		#estimate noise statistics
		if args.locResMap is None: #if no local Resolution map is given,don't do any filtration
			
			mean, var, _ = estimateNoiseFromMap(mapData, wn, boxCoord);
						
			#if varianceMap is given, use it
			if args.varianceMap is not None:
				varMap = mrcfile.open(args.varianceMap, mode='r+');
				var = np.copy(varMap.data);
			
			#if meanMap is given, use it
			if args.meanMap is not None:
				meanMap = mrcfile.open(args.meanMap, mode='r+');
				mean = np.copy(meanMap.data);

			if np.isscalar(mean) and np.isscalar(var):	
				output = "Estimated noise statistics: mean: " + repr(mean) + " and variance: " + repr(var); 
			else:
				output = "Using user provided noise statistics"; 

			print(output);

		else: #do localFiltration and estimate statistics from this map
			locResMap = mrcfile.open(args.locResMap, mode='r+');
			locResMapData = np.copy(locResMap.data);
			mapData, mean, var, ECDF = localFiltration(mapData, locResMapData, apix, True, wn, boxCoord, args.mask, maskData, ECDF);		
		
			locFiltMap = mrcfile.new(splitFilename[0] + '_locFilt.mrc', overwrite=True)
			mapData = np.float32(mapData); 
			locFiltMap.set_data(mapData);
			locFiltMap.voxel_size = apix;
			locFiltMap.close();	

		#calculate the qMap
		qMap = calcQMap(mapData, mean, var, ECDF, wn, boxCoord, circularMaskData, method, testProc);	

		#if a explicit thresholding is wished, do so
		if args.fdr is not None:
			
			fdr = args.fdr;
		
			#threshold the qMap
			binMap = binarizeMap(qMap, fdr);
		
			#apply the thresholded qMap to data
			maskedMap = np.multiply(binMap, mapData);
			minMapValue = np.min(maskedMap[np.nonzero(maskedMap)]);

			maskedMap = np.multiply(maskedMap, maskData);

			if args.locResMap is None: #if no local Resolution map is give, then give the correspoding threshold, not usefule with local filtration 
				output = "Calculated map threshold: " + repr(minMapValue) + " at a FDR of " + repr(fdr);
				print(output);	   
	
			binMapMRC = mrcfile.new(splitFilename[0] + '_FDR' + str(fdr) + '_binMap.mrc', overwrite=True);
			binMap = np.float32(binMap);
			binMapMRC.set_data(binMap);
			binMapMRC.voxel_size = apix;
			binMapMRC.close();

			maskedMapMRC = mrcfile.new(splitFilename[0] + '_FDR'+ str(fdr) + '_maskedMap.mrc', overwrite=True);
			maskedMap = np.float32(maskedMap);
			maskedMapMRC.set_data(maskedMap);
			maskedMapMRC.voxel_size = apix;
			maskedMapMRC.close();
		
		else:
			#threshold the qMap
                        fdr = 0.01;
			binMap = binarizeMap(qMap, fdr);

                        #apply the thresholded qMap to data
                        maskedMap = np.multiply(binMap, np.copy(mapData));
                        minMapValue = np.min(maskedMap[np.nonzero(maskedMap)]);

                        if args.locResMap is None: #if no local Resolution map is give, then give the correspoding threshold, not usefule with local filtration
                        	output = "Calculated map threshold: " + repr(minMapValue) + " at a FDR of " + repr(fdr);
                        	print(output);
	
		#invert qMap for visualization tools
		confidenceMap = np.subtract(1.0, qMap);
	
		#apply lowpass-filtered mask to maps
		confidenceMap = np.multiply(confidenceMap, circularMaskData);
			
		#write the confidence Maps
		confidenceMapMRC = mrcfile.new(splitFilename[0] + '_confidenceMap.mrc', overwrite=True);
		confidenceMap = np.float32(confidenceMap);
		confidenceMapMRC.set_data(confidenceMap);
		confidenceMapMRC.voxel_size = apix;
		confidenceMapMRC.close();

		end = time.time();
		totalRuntime = end -start;

		printSummary(args, totalRuntime);

if (__name__ == "__main__"):
	main()
	







