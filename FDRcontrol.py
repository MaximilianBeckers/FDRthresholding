#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group (2017)

#import some stuff
from confidenceMapUtil import mapUtil, locscaleUtil, FDRutil
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
cmdl_parser.add_argument('-p', '--apix', metavar="apix",  type=float, required=False,
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

		locscaleUtil.launch_amplitude_scaling(args);
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
				map1 = mrcfile.open(args.em_map, mode='r');
				args.apix = map1.voxel_size.x;
				halfMapData1 = np.copy(map1.data);
				
				map2 = mrcfile.open(args.halfmap2, mode='r');
				halfMapData2 = np.copy(map2.data);

				mapData = (halfMapData1 + halfMapData2)*0.5;
				halfMapData1 = 0;
				halfMapData2 = 0;
                                
		else:
			#load single map
			filename = args.em_map;
			map = mrcfile.open(filename, mode='r');
			args.apix = float(map.voxel_size.x);
			mapData = np.copy(map.data);
			
			
		#set output filename
		if args.outputFilename is not None:
			splitFilename = os.path.splitext(os.path.basename(args.outputFilename));
		else:
			splitFilename = os.path.splitext(os.path.basename(filename));


		#if mask is provided, take it
		if args.mask is not None:	
			mask = mrcfile.open(args.mask, mode='r');
			maskData = np.copy(mask.data);
		else:
			maskData = None;

		#if varianceMap is given, use it
		if args.varianceMap is not None:
			varMap = mrcfile.open(args.varianceMap, mode='r');
			varMapData = np.copy(varMap.data);
		else:
			varMapData = None;
			
		#if meanMap is given, use it
		if args.meanMap is not None:
			meanMap = mrcfile.open(args.meanMap, mode='r');
			meanMapData = np.copy(meanMap.data);
		else:
			meanMapData = None;

		#if local resolutions are given, use them
		if args.locResMap is not None:
			locResMap = mrcfile.open(args.locResMap, mode='r');
			locResMapData = np.copy(locResMap.data);
		else:
			locResMapData = None;

		#run the actual analysis
		confidenceMap, locFiltMap, binMap, maskedMap = FDRutil.calculateConfidenceMap(mapData, args.apix, args.noiseBox, args.testProc, args.ecdf, args.lowPassFilter, args.method, args.window_size, locResMapData, meanMapData, varMapData, args.fdr);
			
		if locFiltMap is not None:
			locFiltMapMRC = mrcfile.new(splitFilename[0] + '_locFilt.mrc', overwrite=True);
			locFiltMap = np.float32(locFiltMap);
			locFiltMapMRC.set_data(locFiltMap);
			locFiltMapMRC.voxel_size = args.apix;
			locFiltMapMRC.close();

		if binMap is not None:
			binMapMRC = mrcfile.new(splitFilename[0] + '_FDR' + str(args.fdr) + '_binMap.mrc', overwrite=True);
			binMap = np.float32(binMap);
			binMapMRC.set_data(binMap);
			binMapMRC.voxel_size = args.apix;
			binMapMRC.close();

		if maskedMap is not None:
			maskedMapMRC = mrcfile.new(splitFilename[0] + '_FDR'+ str(args.fdr) + '_maskedMap.mrc', overwrite=True);
			maskedMap = np.float32(maskedMap);
			maskedMapMRC.set_data(maskedMap);
			maskedMapMRC.voxel_size = args.apix;
			maskedMapMRC.close();
	
			
		#write the confidence Maps
		confidenceMapMRC = mrcfile.new(splitFilename[0] + '_confidenceMap.mrc', overwrite=True);
		confidenceMap = np.float32(confidenceMap);
		confidenceMapMRC.set_data(confidenceMap);
		confidenceMapMRC.voxel_size = args.apix;
		confidenceMapMRC.close();

		end = time.time();
		totalRuntime = end -start;

		mapUtil.printSummary(args, totalRuntime);

if (__name__ == "__main__"):
	main()
	







