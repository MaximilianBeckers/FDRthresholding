from EMAN2 import EMData, EMNumPy, Util, XYData, Region
import numpy as np
import argparse, math, os, sys
from FDRutil import *
from mapUtil import *
import time
     
def setup_test_data(voldim=30, size=10):
    from sparx import model_gauss
    emmap = model_gauss(size, voldim, voldim, voldim)
    modmap = EMData()
    modmap.set_size(voldim, voldim, voldim)
    modmap.process_inplace("testimage.noise.gauss", {"sigma":1, "seed":99})
    mask = model_square(size, voldim, voldim, voldim)
    
    return emmap, modmap, mask

def setup_test_data_to_files(emmap_name='emmap.mrc', modmap_name='modmap.mrc', mask_name='mask.mrc'):
    """
    >>> emmap_name, modmap_name, mask_name = setup_test_data_to_files()
    >>> import subprocess
    >>> n = subprocess.call(simple_cmd.split())
    >>> scaled_vol = get_image('scaled.mrc')
    >>> np.copy(EMNumPy.em2numpy(scaled_vol))[scaled_vol.get_xsize() / 2][scaled_vol.get_ysize() / 2]
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.12524424,  0.15562208,  0.18547297,  0.24380369,  0.31203741,
            0.46546721,  0.47914436,  0.31334871,  0.28510684,  0.21345402,
            0.17892323,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ], dtype=float32)
    >>> n = [os.remove(each_file) for each_file in [emmap_name, modmap_name, mask_name, 'scaled.mrc']]
    """
    emmap, modmap, mask = setup_test_data()

    emmap.write_image(emmap_name)
    modmap.write_image(modmap_name)
    mask.write_image(mask_name)
    
    return emmap_name, modmap_name, mask_name

def compute_radial_amplitude_distribution(map, apix):
    data = map.do_fft()
    radial_average = data.calc_radial_dist(map.get_xsize() / 2, 0, 1.0, 0)
    return radial_average
  
def set_zero_origin_and_pixel_size(map, apix):
    map['MRC.nxstart'] = 0
    map['MRC.nystart'] = 0
    map['MRC.nzstart'] = 0
    map.set_attr("apix_x", apix)
    map.set_attr("apix_y", apix)
    map.set_attr("apix_z", apix)
    return map
  
def set_radial_amplitude_distribution(map, amplitude_distribution, apix):
   
    data = map.do_fft()
    frequency_range = np.arange(0, (1 / (2 * apix)), (1.0 / (apix * map.get_xsize())))

    frequency_range = np.ndarray.tolist(frequency_range[0:len(amplitude_distribution)])
	
    sf = XYData()
    sf.set_xy_list(frequency_range, amplitude_distribution)
    data.process_inplace("filter.setstrucfac", {"apix":apix, "strucfac":sf})
      
    mapSharp = data.do_ift()

    return mapSharp
 
def calculateSharpeningFactors(mapObs, mapSharp, apix):
    
    #calculate the two fourier amplitude profiles
    dataObs = mapObs.do_fft()
    radial_averageObs = np.asarray(dataObs.calc_radial_dist(mapObs.get_xsize() / 2, 0, 1.0, 0))

    dataSharp = mapSharp.do_fft()
    radial_averageSharp = np.asarray(dataSharp.calc_radial_dist(mapSharp.get_xsize() / 2, 0, 1.0, 0))

    np.seterr(divide='ignore', invalid='ignore') #no error for division by zero for the following division
    
    #calculate the sharpening factors by division
    sharpFac = np.divide(radial_averageSharp, radial_averageObs)
    sharpFac[ ~ np.isfinite( sharpFac )] = 0 #handle division by zero

    return sharpFac

def get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn):
    mask = np.copy(EMNumPy.em2numpy(mask))
    nk, nj, ni = mask.shape

    kk, jj, ii = np.indices((mask.shape))
    kk_flat = kk.ravel()
    jj_flat = jj.ravel()
    ii_flat = ii.ravel()
    
    mask_bin = np.array(mask.ravel(), dtype=np.bool)
    indices = np.arange(mask.size)
    masked_indices = indices[mask_bin]
    cropped_indices = indices[(wn / 2 <= kk_flat) & (kk_flat < (nk - wn / 2)) &
                              (wn / 2 <= jj_flat) & (jj_flat < (nj - wn / 2)) &
                              (wn / 2 <= ii_flat) & (ii_flat < (ni - wn / 2))]
                                     
    cropp_n_mask_ind = np.intersect1d(masked_indices, cropped_indices)
    
    xyz_locs = np.column_stack((kk_flat[cropp_n_mask_ind], jj_flat[cropp_n_mask_ind], ii_flat[cropp_n_mask_ind]))
    
    return xyz_locs, cropp_n_mask_ind, mask.shape

def get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, locFilt, locResMap, boxCoord):
    sharpened_vals = np.array([], dtype=np.float32)
    sharpened_mean_vals = np.array([], dtype=np.float32)
    sharpened_var_vals = np.array([], dtype=np.float32)
    
    reg = Region(wn, wn, wn, wn, wn, wn)

    central_pix = int(round(wn / 2.0))
    sizeMap = np.array([emmap.get_xsize(), emmap.get_ysize(), emmap.get_zsize()])
    center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]])

    emmapData = np.copy(EMNumPy.em2numpy(emmap))
    
    if boxCoord == 0:
        noiseMapData = emmapData[int(center[0]-0.5*wn):(int(center[0]-0.5*wn) + wn), int(0.02*wn+wn):(int(0.02*wn+wn) + wn), (int(center[2]-0.5*wn)):(int((center[2]-0.5*wn) + wn))]
    else:
        noiseMapData = emmapData[int(boxCoord[0]-0.5*wn+wn):(int(boxCoord[0]-0.5*wn+wn) + wn), int(boxCoord[1]+wn-0.5*wn):(int(boxCoord[1]+wn-0.5*wn) + wn), (int(boxCoord[2]+wn-0.5*wn)):(int((boxCoord[2]+wn-0.5*wn)+wn))]

    noiseMap = EMNumPy.numpy2em(np.copy(noiseMapData))
    noise_radial_average = compute_radial_amplitude_distribution(noiseMap, apix)
    noise_radial_average_np = np.asarray(noise_radial_average)
	
	#######
    nx, ny, nz = noiseMapData.shape;
    noise_image_sample = np.zeros(( nx, ny)); 
    counterNoiseSample = 0;


    for k, j, i in (masked_xyz_locs - wn / 2):
		
        reg = Region(i, j, k, wn, wn, wn)

        emmap_wn = emmap.get_clip(reg)
        modmap_wn = modmap.get_clip(reg)

        emmap_radial_average = compute_radial_amplitude_distribution(emmap_wn, apix)
        mod_radial_average = compute_radial_amplitude_distribution(modmap_wn, apix)
		
        map_b_sharpened = set_radial_amplitude_distribution(emmap_wn, mod_radial_average, apix)
      	
        #get scaling factors		
        sharpFactors = calculateSharpeningFactors(emmap_wn, map_b_sharpened, apix)		

        #scale noise window with the same scaling factors
        noiseMapFFT = noiseMap.do_fft()
        noiseMapFFT.apply_radial_func(0.0, (1.0 / (apix * noiseMap.get_xsize())), sharpFactors.tolist())
        mapNoise_sharpened = noiseMapFFT.do_ift()
        
        	
        if locFilt == True:
            tmpRes = round(apix/locResMap[k, j, i],3)
            				
            mapNoise_sharpened = mapNoise_sharpened.process("filter.lowpass.tanh", {"cutoff_abs": tmpRes, "fall_off": 0.1})
            map_b_sharpened = map_b_sharpened.process("filter.lowpass.tanh", {"cutoff_abs": tmpRes, "fall_off": 0.1})     
			 
            #calculate noise statistics	  
            map_noise_sharpened_data = np.copy(EMNumPy.em2numpy(mapNoise_sharpened))	   
            mean = np.mean(map_noise_sharpened_data)
            var = np.var(map_noise_sharpened_data)
            if tmpRes == round(apix/100.0, 3):  
                mean = 0.0
                var = 0.0
        else:
            #calculate noise statistics
            map_noise_sharpened_data = np.copy(EMNumPy.em2numpy(mapNoise_sharpened))
            mean = np.mean(map_noise_sharpened_data)
            var = np.var(map_noise_sharpened_data)


        #append values to the sharpened values
        sharpened_vals = np.append(sharpened_vals, map_b_sharpened[central_pix, central_pix, central_pix])
        sharpened_mean_vals = np.append(sharpened_mean_vals, mean)
        sharpened_var_vals = np.append(sharpened_var_vals, var)     
        
    return sharpened_vals, sharpened_mean_vals, sharpened_var_vals

def put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape):
    map_scaled = np.zeros(np.prod(map_shape))
    map_scaled[masked_indices] = sharpened_vals
    map_scaled = map_scaled.reshape(map_shape)
    
    map_scaled = EMNumPy.numpy2em(np.copy(map_scaled))
    
    return map_scaled

def run_window_function_including_scaling(emmap, modmap, mask, wn, apix, locFilt, locResMap, boxCoord):
    """
    >>> emmap, modmap, mask = setup_test_data()
    >>> scaled_vol = run_window_function_including_scaling(emmap,modmap,mask,wn=10,apix=1.0)
    >>> np.copy(EMNumPy.em2numpy(scaled_vol))[scaled_vol.get_xsize() / 2][scaled_vol.get_ysize() / 2]
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.12524424,  0.15562208,  0.18547297,  0.24380369,  0.31203741,
            0.46546721,  0.47914436,  0.31334871,  0.28510684,  0.21345402,
            0.17892323,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ], dtype=float32)
    """
    masked_xyz_locs, masked_indices, map_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)
 
    sharpened_vals, sharpened_mean_vals, sharpened_var_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, locFilt, locResMap, boxCoord)
     
    map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)
    mean_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)	 
    var_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)

    return map_scaled, mean_map_scaled, var_map_scaled 
   
def split_sequence_evenly(seq, size):
    """
    >>> split_sequence_evenly(list(range(9)), 4)
    [[0, 1], [2, 3, 4], [5, 6], [7, 8]]
    >>> split_sequence_evenly(list(range(18)), 4)
    [[0, 1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12, 13], [14, 15, 16, 17]]
    """
    newseq = []
    splitsize = 1.0 / size * len(seq)
    for i in range(size):
        newseq.append(seq[int(round(i * splitsize)):int(round((i + 1) * splitsize))])
    return newseq
    
def merge_sequence_of_sequences(seq):
    """
    >>> merge_sequence_of_sequences([list(range(9)), list(range(3))])
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2]
    >>> merge_sequence_of_sequences([list(range(9)), [], list(range(3))])
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2]
    """
    newseq = [number for sequence in seq for number in sequence]
    
    return newseq
    
def compute_padding_average(map, mask):
    volume_stats_outside_mask = Util.infomask(map, mask, False)
    average_padding_intensity = volume_stats_outside_mask[0]
    
    return average_padding_intensity
    
def check_for_window_bleeding(mask,wn):
    masked_xyz_locs, masked_indices, mask_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, 0)
    
    zs, ys, xs = masked_xyz_locs.T
    nk, nj, ni = mask_shape

    if xs.min() < wn / 2 or xs.max() > (ni - wn / 2) or \
    ys.min() < wn / 2 or ys.max() > (nj - wn / 2) or \
    zs.min() < wn / 2 or zs.max() > (nk - wn / 2):
        window_bleed = True
    else:
        window_bleed = False
    
    return window_bleed 	
    
def prepare_mask_and_maps_for_scaling(args):
    
    from sparx import get_image, binarize, model_square	
    emmap = get_image(args.em_map)
    modmap = get_image(args.model_map)

    if args.mask is None:
        mask = EMData()
        xsize, ysize, zsize = emmap.get_xsize(), emmap.get_ysize(), emmap.get_zsize()
        mask.set_size(xsize, ysize, zsize)
        mask.to_zero()
        if xsize == ysize and xsize == zsize and ysize == zsize:
            sphere_radius = xsize // 2
            mask.process_inplace("testimage.circlesphere", {"radius":sphere_radius})
        else:
            mask += 1
            mask = Util.window(mask, xsize - 1, ysize - 1, zsize -1)
            mask = Util.pad(mask, xsize, ysize, zsize, 0, 0, 0, '0')
    elif args.mask is not None:
        mask = binarize(get_image(args.mask), 0.5)
        
    if args.window_size is None:
        wn = int(round(7 * 3 * args.apix)) # set default window size to 7 times average resolution
    elif args.window_size is not None:
        wn = int(math.ceil(args.window_size / 2.) * 2)
    
    if args.FDRmethod is not None:
        FDRmethod = args.FDRmethod
    else:
        FDRmethod = 'BY'
	
    if args.noiseBox is not None:
        boxCoord = args.noiseBox;
    else:
        boxCoord = 0;

    if args.locResMap is not None:
        locResMap = get_image(args.locResMap)

    window_bleed_and_pad = check_for_window_bleeding(mask, wn)
    if window_bleed_and_pad:
        pad_int_emmap = compute_padding_average(emmap, mask)
        pad_int_modmap = compute_padding_average(modmap, mask)

        map_shape = [(emmap.get_xsize() + wn), (emmap.get_ysize() + wn), (emmap.get_zsize() + wn)] 
        emmap = Util.pad(emmap, map_shape[0], map_shape[1], map_shape[2], 0, 0, 0, 'pad_int_emmap')
        modmap = Util.pad(modmap, map_shape[0], map_shape[1], map_shape[2], 0, 0, 0, 'pad_int_modmap')
        mask = Util.pad(mask, map_shape[0], map_shape[1], map_shape[2], 0, 0, 0, '0')
        if args.locResMap is not None:
            locResMap = Util.pad(locResMap, map_shape[0], map_shape[1], map_shape[2], 0, 0, 0, '100.0');

    #if wished so, do local filtration                      
    if args.locResMap is not None:
        locResMapData = np.copy(EMNumPy.em2numpy(locResMap))  
        locResMapData[locResMapData == 0.0] = 100.0;
        locFilt = True    
    else:
        locFilt = False 
        locResMapData = np.ones([xsize, ysize, zsize]);

    return emmap, modmap, mask, wn, window_bleed_and_pad, FDRmethod, locFilt, locResMapData, boxCoord

def run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, apix, locFilt, locResMap, boxCoord):
    """
    >>> emmap_name, modmap_name, mask_name = setup_test_data_to_files()
    >>> import subprocess
    >>> n = subprocess.call(mpi_cmd.split())
    >>> scaled_vol = get_image('scaled.mrc')
    >>> np.copy(EMNumPy.em2numpy(scaled_vol))[scaled_vol.get_xsize() / 2][scaled_vol.get_ysize() / 2]
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.12524424,  0.15562208,  0.18547297,  0.24380369,  0.31203741,
            0.46546721,  0.47914436,  0.31334871,  0.28510684,  0.21345402,
            0.17892323,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ], dtype=float32)
    >>> n = [os.remove(each_file) for each_file in [emmap_name, modmap_name, mask_name, 'scaled.mrc']]
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    if rank == 0:
        
        print("*****************************")
        print("********* LocScale **********")
        print("*****************************")


        masked_xyz_locs, masked_indices, map_shape = \
        get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)
         
        zs, ys, xs = masked_xyz_locs.T
        zs = split_sequence_evenly(zs, size)
        ys = split_sequence_evenly(ys, size)
        xs = split_sequence_evenly(xs, size)
    else:
        zs = None
        ys = None
        xs = None

    comm.barrier()

    zs = comm.scatter(zs, root=0)
    ys = comm.scatter(ys, root=0)
    xs = comm.scatter(xs, root=0)
 
    masked_xyz_locs = np.column_stack((zs, ys, xs))
 
    sharpened_vals, sharpened_mean_vals, sharpened_var_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, locFilt, locResMap, boxCoord)

    comm.barrier()	
    sharpened_vals = comm.gather(sharpened_vals, root=0)
    sharpened_mean_vals = comm.gather(sharpened_mean_vals, root=0)
    sharpened_var_vals = comm.gather(sharpened_var_vals, root=0)		
 
    if rank == 0:
        sharpened_vals = merge_sequence_of_sequences(sharpened_vals)
        sharpened_mean_vals = merge_sequence_of_sequences(sharpened_mean_vals)
        sharpened_var_vals = merge_sequence_of_sequences(sharpened_var_vals)

        map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_vals),
        masked_indices, map_shape)

        mean_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_mean_vals),
        masked_indices, map_shape)

        var_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_var_vals),
        masked_indices, map_shape)
    else:
        map_scaled = None
        mean_map_scaled = None
        var_map_scaled = None		
 
    return map_scaled, mean_map_scaled, var_map_scaled, rank
  
def write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol, filename):
    LocScaleVol = set_zero_origin_and_pixel_size(LocScaleVol, args.apix)
    if window_bleed_and_pad:
        map_shape = [(LocScaleVol.get_xsize() - wn), (LocScaleVol.get_ysize() - wn), (LocScaleVol.get_zsize() - wn)]
        LocScaleVol = Util.window(LocScaleVol, map_shape[0], map_shape[1], map_shape[2])

    LocScaleVol.write_image(filename)

    return LocScaleVol

def launch_amplitude_scaling(args):    	


    startTime = time.time()
    emmap, modmap, mask, wn, window_bleed_and_pad, FDRmethod, locFilt, locResMap, boxCoord = prepare_mask_and_maps_for_scaling(args)
     
    meanNoise, varNoise, sample = estimateNoiseFromMap(EMNumPy.em2numpy(emmap), wn, boxCoord)
   
    #set output filenames
    if args.outputFilename is not None:
        splitFilename = os.path.splitext(os.path.basename(args.outputFilename))
    else:	
        splitFilename = os.path.splitext(os.path.basename(args.em_map))
	

    if not args.mpi:

        LocScaleVol, meanVol, varVol = run_window_function_including_scaling(emmap, modmap, mask, wn, args.apix, locFilt, locResMap, boxCoord)
        print("Local amplitude scaling finished.")

        #calculate qMap
        print("Start significance analysis.")
        meanVolData = EMNumPy.em2numpy(meanVol)
        varVolData = EMNumPy.em2numpy(varVol)
        LocScaleVolData = EMNumPy.em2numpy(LocScaleVol)
        maskData = EMNumPy.em2numpy(mask)

        qVolData = calcQMap(LocScaleVolData, meanVolData, varVolData, maskData, FDRmethod, 'rightSided')
        qVolData = np.subtract(np.ones(qVolData.shape), qVolData)
        qVol = EMNumPy.numpy2em(qVolData)

        #write the volumes    
        LocScaleVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol, splitFilename[0] + '_scaled.mrc')
        qVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, qVol, splitFilename[0] + '_locscale' + '_confidenceMap.mrc')
        #meanVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, meanVol, 'meanVol.mrc')
        #varVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, varVol, 'varVol.mrc') 

		#make diagnostic plot and print summary 
        makeDiagnosticPlot(emmap, wn, wn, True, boxCoord);   

        endTime = time.time()
        runTime = endTime - startTime
        printSummary(args, runTime)

    elif args.mpi:


        LocScaleVol, meanVol, varVol, rank = run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, args.apix, locFilt, locResMap, boxCoord)
	
        if rank == 0:
            print("Local amplitude scaling finished.")			

            #calculate qMap
            print("Start significance analysis.")
            meanVolData = EMNumPy.em2numpy(meanVol)
            varVolData = EMNumPy.em2numpy(varVol)
            LocScaleVolData = EMNumPy.em2numpy(LocScaleVol)
            maskData = EMNumPy.em2numpy(mask)

            qVolData = calcQMap(LocScaleVolData, meanVolData, varVolData, maskData, FDRmethod, 'rightSided')
            #invert qMap for visualization tools
            qVolData = np.subtract(np.ones(qVolData.shape), qVolData)
            qVol = EMNumPy.numpy2em(qVolData)

			#write the volumes
            LocScaleVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol, splitFilename[0] + '_scaled.mrc')
            qVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, qVol, splitFilename[0] + '_confidenceMap.mrc')
            meanVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, meanVol, 'meanVol.mrc')
            varVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, varVol, 'varVol.mrc') 

			#make diagnostic plot and print summary 
            makeDiagnosticPlot(emmap, wn, wn, True, boxCoord);   
            
            endTime = time.time()
            runTime = endTime - startTime			
            printSummary(args, runTime)

