#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group (2017)

#import some stuff
from EMAN2 import EMData, EMNumPy
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

#*************************************************************
#****************** Commandline input ************************
#*************************************************************

cmdl_parser = argparse.ArgumentParser(
               prog=sys.argv[0],
               description='*** Analyse density ***',
               formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=30), add_help=True);

cmdl_parser.add_argument('-em', '--em_map', metavar="em_map.mrc",  type=str, required=True,
                         help='Input filename EM map');
cmdl_parser.add_argument('-mask', '--mask',metavar="mask",  type=str,
                         help='Input filename mask');
cmdl_parser.add_argument('-p', '--apix', metavar="apix",  type=float, required=True,
                         help='pixel Size of input map');
cmdl_parser.add_argument('-fdr', '--fdr', metavar="fdr",  type=float, required=False,
                         help='False Discovery Rate');
cmdl_parser.add_argument('-locResMap', metavar="locResMap.mrc", type=str, required=False,
						help='Input local Resolution Map');
cmdl_parser.add_argument('-FDRmethod', metavar="FDRmethod", type=str, required=False,
						help="Input FDR method. 'BY' for Benjamini-Yekutieli or 'BH' for Benjamini-Hochberg");
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

#************************************************************
#********************** main function ***********************
#************************************************************

def main():

	start = time.time();

	#get command line input
	args = cmdl_parser.parse_args();

	#if a model map is given, local amplitude scaling will be done
	if args.model_map is not None:
		launch_amplitude_scaling(args);

	else:	
		#no ampltidue scaling will be done
		print('************************************************');
		print('******* Significance analysis of EM-Maps *******');
		print('************************************************');
	
		filename = args.em_map;
		map = EMData();
		map.read_image(filename);
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

		sizeMap = np.array([map.get_xsize(), map.get_ysize(), map.get_zsize()]);

		if args.lowPassFilter is not None:
			providedRes = apix/float(args.lowPassFilter);
			map.process_inplace("filter.lowpass.tanh", {"cutoff_abs": providedRes, "fall_off": 0.1});

		#handle FDR correction procedure
		if args.FDRmethod is not None:
			FDRmethod = args.FDRmethod;
		else:
			#default is Benjamini-Yekutieli
			FDRmethod = 'BY';

		if args.window_size is not None:
			wn = args.window_size;
			wn = int(wn);
			if wn < 10:
				print("Provided window size is very small. Please think about potential inaccuracies of your noise estimates!");
		else: 
			wn = max(int(0.05*sizeMap[0]),10);

		#generate a circular Mask
		circularMask = EMData();
		circularMask.set_size(sizeMap[0], sizeMap[1], sizeMap[2]);
		circularMask.to_zero();	
		sphere_radius = (np.min(sizeMap) // 2) + 10;
		circularMask.process_inplace("testimage.circlesphere", {"radius":sphere_radius});
		#circularMask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1});			
		circularMaskData = np.copy(EMNumPy.em2numpy(circularMask));
		
		#if mask is provided, take it
		if args.mask is not None:	
			mask = EMData();
			mask.read_image(args.mask);
			maskData = np.copy(EMNumPy.em2numpy(mask));
		else:
			maskData = circularMaskData;
			

		#estimate noise statistics
		if args.locResMap is None: #if no local Resolution map is given,don't do any filtration
			mapData = EMNumPy.em2numpy(map);
			
			if args.mask is None:
				mean, var, _ = estimateNoiseFromMap(mapData, wn, boxCoord);
			else:
				mean, var, _ = estimateNoiseFromMapInsideMask(mapData, maskData);
		
			#if varianceMap is given, use it
			if args.varianceMap is not None:
				varMap = EMData();
				varMap.read_image(args.varianceMap);
				var = np.copy(EMNumPy.em2numpy(varMap));
			
			#if meanMap is given, use it
			if args.meanMap is not None:
				meanMap = EMData();
				meanMap.read_image(args.meanMap);
				mean = np.copy(EMNumPy.em2numpy(meanMap));

			if np.isscalar(mean) and np.isscalar(var):	
				output = "Estimated noise statistics: mean: " + repr(mean) + " and variance: " + repr(var); 
			else:
				output = "Using user provided noise statistics"; 

			print(output);

		else: #do localFiltration and estimate statistics from this map
			locResMap = EMData();
			locResMap.read_image(args.locResMap);
			mapData, mean, var, ECDF = localFiltration(map, locResMap, apix, True, wn, boxCoord, args.mask, maskData, ECDF);		

			locFiltMapEMAN = EMNumPy.numpy2em(mapData); 
			locFiltMapEMAN =  set_zero_origin_and_pixel_size(locFiltMapEMAN, apix);
			locFiltMapEMAN.write_image(splitFilename[0] + '_locFilt.mrc');

		makeDiagnosticPlot(map, wn, 0, False, boxCoord);

		#calculate the qMap
		qMap = calcQMap(mapData, mean, var, ECDF, wn, boxCoord, circularMaskData, FDRmethod, testProc);	

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
	
			binMapEMAN = EMNumPy.numpy2em(binMap);
			binMapEMAN = set_zero_origin_and_pixel_size(binMapEMAN, apix);
			binMapEMAN.write_image(splitFilename[0] + '_FDR' + str(fdr) + '_binMap.mrc');
     
			maskedMapEMAN = EMNumPy.numpy2em(maskedMap);
			maskedMapEMAN = set_zero_origin_and_pixel_size(maskedMapEMAN, apix);
			maskedMapEMAN.write_image(splitFilename[0] + '_FDR'+ str(fdr) + '_maskedMap.mrc');

	
		#invert qMap for visualization tools
		confidenceMap = np.subtract(1.0, qMap);
	
		#apply lowpass-filtered mask to maps
		confidenceMap = np.multiply(confidenceMap, circularMaskData);
			
		#write the confidence Maps
		confidenceMapEMAN = EMNumPy.numpy2em(confidenceMap);
		confidenceMapEMAN = set_zero_origin_and_pixel_size(confidenceMapEMAN, apix);
		confidenceMapEMAN.write_image(splitFilename[0] + '_confidenceMap.mrc');	
	
		end = time.time();
		totalRuntime = end -start;

		printSummary(args, totalRuntime);

if (__name__ == "__main__"):
	main()
	







