#!/programs/x86_64-linux/eman2/2.2/bin/python
#Author: Maximilian Beckers, EMBL Heidelberg, Sachse Group (2017)

#import some stuff
from EMAN2 import *
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
               formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=30), add_help=False)

cmdl_parser.add_argument('-em', '--em_map', metavar="em_map.mrc",  type=str, required=True,
                         help='Input filename EM map')
cmdl_parser.add_argument('-mask', '--mask',metavar="mask",  type=str,
                         help='Input filename mask')
cmdl_parser.add_argument('-p', '--apix', metavar="apix",  type=float, required=True,
                         help='pixel Size of input map')
cmdl_parser.add_argument('-fdr', '--fdr', metavar="fdr",  type=float, required=False,
                         help='False Discovery Rate')
cmdl_parser.add_argument('-locResMap', metavar="locResMap.mrc", type=str, required=False,
						help='Input local Resolution Map')
cmdl_parser.add_argument('-FDRmethod', metavar="FDRmethod", type=str, required=False,
						help="Input FDR method. 'BY' for Benjamini-Yekutieli or 'BH' for Benjamini-Hochberg")
cmdl_parser.add_argument('-mm', '--model_map', metavar="model_map.mrc", type=str, required=False,
						help="Input model map for model based amplitude scaling")
cmdl_parser.add_argument('-w', '--window_size', metavar="windowSize", type=float, required=False,
						help="Input window size for local Amplitude scaling and background noise estimation")
cmdl_parser.add_argument('-mpi', action='store_true', default=False,
						help="Set this flag if MPI should be used for the local amplitude scaling")
cmdl_parser.add_argument('-o', '--output_filename', metavar="output.mrc", type=str, required=False,
						help="Name of the output")


#************************************************************
#********************** main function ***********************
#************************************************************

def main():

	start = time.time()

	#get command line input
	args = cmdl_parser.parse_args()

	#if a model map is given, local amplitude scaling will be done
	if args.model_map is not None:
		launch_amplitude_scaling(args)

	else:	
		#no ampltidue scaling will be done
		print('************************************************')
		print('******* Significance analysis of EM-Maps *******')
		print('************************************************')
	
		filename = args.em_map
		map = EMData()
		map.read_image(filename)
		apix = args.apix

		#set the output filenames
		if args.output_filename is not None:
			splitFilename = os.path.splitext(os.path.basename(args.output_filename))
		else:
			splitFilename = os.path.splitext(os.path.basename(filename))

		sizeMap = np.array([map.get_xsize(), map.get_ysize(), map.get_zsize()])

		#handle FDR correction procedure
		if args.FDRmethod is not None:
			FDRmethod = args.FDRmethod
		else:
			#default is Benjamini-Yekutieli
			FDRmethod = 'BY'

		if args.window_size is not None:
			wn = args.window_size
		else: 
			wn = int(0.05*sizeMap[0])

		#handle masking options
		if args.mask is None: #if no mask is given, generate a circular one
			mask = EMData()
			mask.set_size(sizeMap[0], sizeMap[1], sizeMap[2])
			mask.to_zero()	
			sphere_radius = np.min(sizeMap) // 2
			mask.process_inplace("testimage.circlesphere", {"radius":sphere_radius})
			mask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})			
			maskData = EMNumPy.em2numpy(mask)

		else: #else use the provided mask
			mask = EMData()
			mask.read_image(args.mask)
			maskData = EMNumPy.em2numpy(mask)

		#estimate noise statistics
		if args.locResMap is None: #if no local Resolution map is given,don't do any filtration
			mapData = EMNumPy.em2numpy(map)
			mean, var = estimateNoiseFromMap(mapData, wn)
		else: #do localFiltration and estimate statistics from this map
			locResMap = EMData()
			locResMap.read_image(args.locResMap)
			mapData, mean, var = localFiltration(map, locResMap, apix, True, wn)		

			locFiltMapEMAN = EMNumPy.numpy2em(mapData) 
			locFiltMapEMAN.write_image(splitFilename[0] + '_locFilt.mrc')

		makeDiagnosticPlot(map, wn, 0, False)

		#calculate the qMap
		qMap = calcQMap(mapData, mean, var, maskData, FDRmethod, 'rightSided')	
	
		if args.fdr is not None:
			
			fdr = args.fdr
		
			#threshold the qMap
			binMap = binarizeMap(qMap, fdr)
		
			#apply the thresholded qMap to data
			maskedMap = np.multiply(binMap, mapData)
			minMapValue = np.min(maskedMap[np.nonzero(maskedMap)])

			binMap = np.multiply(binMap, maskData)
			maskedMap = np.multiply(maskedMap, maskData)

			if args.locResMap is None: #if no local Resolution map is give, then give the correspoding threshold, not usefule with local filtration 
				output = "Calculated map threshold: " + repr(minMapValue) + " at a FDR of " + repr(fdr)
				print(output)	   
	
			binMapEMAN = EMNumPy.numpy2em(binMap)
			binMapEMAN.write_image(splitFilename[0] + '_FDR' + str(fdr) + '_binMap.mrc')
     
			maskedMapEMAN = EMNumPy.numpy2em(maskedMap)
			maskedMapEMAN.write_image(splitFilename[0] + '_FDR'+ str(fdr) + '_maskedMap.mrc')

		#write the qMaps
		#qMapEMAN = EMNumPy.numpy2em(qMap)
		#qMapEMAN.write_image(splitFilename[0] + '_qMap.mrc')
	
		#invert qMap for visualization tools
		qMap = np.subtract(np.ones(sizeMap), qMap)
	
		#apply lowpass-filtered mask to maps
		qMap = np.multiply(qMap, maskData)
	
		#write the qMaps
		qMapEMAN = EMNumPy.numpy2em(qMap)
		qMapEMAN.write_image(splitFilename[0] + 'confidenceMap.mrc')	
	
		end = time.time()
		totalRuntime = end -start

		printSummary(args, totalRuntime)

if (__name__ == "__main__"):
	main()
	







