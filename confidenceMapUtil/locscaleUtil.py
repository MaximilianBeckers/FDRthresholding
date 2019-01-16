import numpy as np
from FDRutil import *
from mapUtil import *
import mrcfile
import argparse, math, os, sys
from argparse import RawTextHelpFormatter
import time
  
def compute_padding_average(vol, mask):
    mask = (mask > 0.5).astype(np.int8)
    average_padding_intensity = np.mean(np.ma.masked_array(vol, mask))
    return average_padding_intensity

def pad_or_crop_volume(vol, dim_pad=None, pad_value = None, crop_volume=False):   
    if (dim_pad == None):
        return vol
    else:
        dim_pad = np.round(np.array(dim_pad)).astype('int')

        if pad_value == None:
            pad_value = 0

        if (dim_pad[0] <= vol.shape[0] or dim_pad[1] <= vol.shape[1] or dim_pad[2] <= vol.shape[2]): 
            crop_volume = True
        
        if crop_volume:
            crop_vol = vol[vol.shape[0]/2-dim_pad[0]/2:vol.shape[0]/2+dim_pad[0]/2+dim_pad[0]%2, :, :]
            crop_vol = crop_vol[:, vol.shape[1]/2-dim_pad[1]/2:vol.shape[1]/2+dim_pad[1]/2+dim_pad[1]%2, :]
            crop_vol = crop_vol[:, :, vol.shape[2]/2-dim_pad[2]/2:vol.shape[2]/2+dim_pad[2]/2+dim_pad[2]%2]
            
            return crop_vol
            
        else:
            pad_vol = np.pad(vol, ((dim_pad[0]/2-vol.shape[0]/2, dim_pad[0]/2-vol.shape[0]/2+dim_pad[0]%2), (0,0), (0,0) ), 'constant', constant_values=(pad_value,))
            pad_vol = np.pad(pad_vol, ((0,0), (dim_pad[1]/2-vol.shape[1]/2, dim_pad[1]/2-vol.shape[1]/2+dim_pad[1]%2 ), (0,0)), 'constant', constant_values=(pad_value,))
            pad_vol = np.pad(pad_vol, ((0,0), (0,0), (dim_pad[2]/2-vol.shape[2]/2, dim_pad[2]/2-vol.shape[2]/2+dim_pad[2]%2)), 'constant', constant_values=(pad_value,))
            
            return pad_vol

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

def get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn):

    mask = np.copy(mask);
    nk, nj, ni = mask.shape;

    kk, jj, ii = np.indices((mask.shape));
    kk_flat = kk.ravel();
    jj_flat = jj.ravel();
    ii_flat = ii.ravel();
    
    mask_bin = np.array(mask.ravel(), dtype=np.bool);
    indices = np.arange(mask.size);
    masked_indices = indices[mask_bin];
    cropped_indices = indices[(wn / 2 <= kk_flat) & (kk_flat < (nk - wn / 2)) &
                              (wn / 2 <= jj_flat) & (jj_flat < (nj - wn / 2)) &
                              (wn / 2 <= ii_flat) & (ii_flat < (ni - wn / 2))];
                                     
    cropp_n_mask_ind = np.intersect1d(masked_indices, cropped_indices);
    
    xyz_locs = np.column_stack((kk_flat[cropp_n_mask_ind], jj_flat[cropp_n_mask_ind], ii_flat[cropp_n_mask_ind]));

    return xyz_locs, cropp_n_mask_ind, mask.shape;

def prepare_mask_and_maps_for_scaling(args):

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

            emmap = (halfMapData1 + halfMapData2)*0.5;
            halfMapData1 = 0;
            halfMapData2 = 0;
                                
    else:
            #load single map
            filename = args.em_map;
            map = mrcfile.open(filename, mode='r+');
            emmap = np.copy(map.data);

    modmap = np.copy(mrcfile.open(args.model_map).data);
    
    if args.mask is None:
        mask = np.zeros(emmap.shape);
        
        if mask.shape[0] == mask.shape[1] and mask.shape[0] == mask.shape[2] and mask.shape[1] == mask.shape[2]:
            rad = (mask.shape[0] // 2) ;
            z,y,x = np.ogrid[-rad: rad+1, -rad: rad+1, -rad: rad+1];
            mask = (x**2+y**2+z**2 <= rad**2).astype(np.int_).astype(np.int8);
            mask = pad_or_crop_volume(mask,emmap.shape);
            mask = (mask > 0.5).astype(np.int8);
        else:
            mask += 1;
            mask = mask[0:mask.shape[0]-1, 0:mask.shape[1]-1, 0:mask.shape[2]-1];
            mask = pad_or_crop_volume(emmap, (emmap.shape), pad_value=0);
		
    elif args.mask is not None:
        mask = (mrcfile.open(args.mask).data > 0.5).astype(np.int8);
                
    if args.window_size_locscale is None:
        wn_locscale = int(round(7 * 3 * args.apix)); # set default window size to 7 times average resolution
    elif args.window_size_locscale is not None:
        wn_locscale = int(math.ceil(args.window_size_locscale / 2.) * 2);
   
    if args.window_size is None:
        wn = wn_locscale;
    elif args.window_size is not None: 
        wn = int(math.ceil(args.window_size / 2.) * 2);

    if args.method is not None:
        method = args.method;
    else:
        method = 'BY';
	
    if args.noiseBox is not None:
        boxCoord = args.noiseBox;
    else:
        boxCoord = 0;

    if args.locResMap is not None:
        locResMap = mrcfile.open(args.locResMap).data;
   
    window_bleed_and_pad = check_for_window_bleeding(mask, wn_locscale);
    if window_bleed_and_pad:
        pad_int_emmap = compute_padding_average(emmap, mask);
        pad_int_modmap = compute_padding_average(modmap, mask);
        map_shape = [(emmap.shape[0] + wn_locscale), (emmap.shape[1] + wn_locscale), (emmap.shape[2] + wn_locscale)];
        emmap = pad_or_crop_volume(emmap, map_shape, pad_int_emmap);
        modmap = pad_or_crop_volume(modmap, map_shape, pad_int_modmap);
        mask = pad_or_crop_volume(mask, map_shape, 0);
        if args.locResMap is not None:
            locResMap = pad_or_crop_volume(locResMap, map_shape, 100.0);
    
    #if wished so, do local filtration                      
    if args.locResMap is not None:
        locResMapData = np.copy(locResMap);  
        locResMapData[locResMapData == 0.0] = 100.0;
        locResMapData[locResMapData >= 100.0] = 100.0;
        locFilt = True;    
    else:
        locFilt = False;
        locResMapData = np.ones(emmap.shape);

    return emmap, modmap, mask, wn, wn_locscale, window_bleed_and_pad, method, locFilt, locResMapData, boxCoord;

def compute_radial_profile(volFFT, frequencyMap):

    dim = volFFT.shape;
    ps = np.abs(volFFT);
    frequencies = np.linspace(0, 0.5, int(math.ceil(dim[0]/2.0)));
    bins = np.digitize(frequencyMap, frequencies);
    bins = bins - 1;
    radial_profile = np.bincount(bins.ravel(), ps.ravel()) / np.bincount(bins.ravel())

    return radial_profile, frequencies;

def compute_scale_factors(em_profile, ref_profile):

    np.seterr(divide='ignore', invalid='ignore'); #no error for division by zero 
    #scale_factor = (ref_profile**2/em_profile**2);
    #scale_factor[ ~ np.isfinite( scale_factor )] = 0 #handle division by zero
    #scale_factor = np.sqrt(scale_factor);
    scale_factor = np.abs(ref_profile/em_profile);    
    scale_factor[ ~ np.isfinite( scale_factor )] = 0; #handle division by zero

    return scale_factor;

def set_radial_profile(volFFT, scaleFactors, frequencies, frequencyMap, shape): 
    
    scalingMap = np.interp(frequencyMap, frequencies, scaleFactors, right=1.0);	
    scaledMapFFT = scalingMap * volFFT;
    scaledMap = np.real(np.fft.irfftn(scaledMapFFT, shape));
    
    return scaledMap, scaledMapFFT;

def calculate_scaled_map(emmap, modmap, mask, wn, wn_locscale, apix, locFilt, locResMap, boxCoord, ecdfBool, stepSize):

    sizeMap = emmap.shape
    sharpened_map = np.zeros(sizeMap);
    sharpened_mean_vals = np.zeros(sizeMap);
    sharpened_var_vals = np.zeros(sizeMap);
    sharpened_ecdf_vals = np.zeros(sizeMap);
    central_pix = int(round(wn_locscale / 2.0));
    center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]]);
    
    #get the background noise sample
    if boxCoord == 0:
        noiseMap = emmap[int(center[0]-0.5*wn):(int(center[0]-0.5*wn) + wn), int(0.02*wn+wn_locscale):(int(0.02*wn+wn_locscale) + wn), (int(center[2]-0.5*wn)):(int((center[2]-0.5*wn) + wn))];
    else:
        noiseMap = emmap[int(boxCoord[0]-0.5*wn +wn_locscale):(int(boxCoord[0]-0.5*wn + wn_locscale) + wn), int(boxCoord[1]-0.5*wn+ wn_locscale):(int(boxCoord[1]-0.5*wn + wn_locscale) + wn), (int(boxCoord[2]-0.5*wn + wn_locscale)):(int((boxCoord[2]-0.5*wn + wn_locscale)+wn))];

    #prepare noise map for scaling
    frequencyMap_noise = calculate_frequency_map(noiseMap);
    noiseMapFFT = np.fft.rfftn(noiseMap);
    noise_profile, frequencies_noise = compute_radial_profile(noiseMapFFT, frequencyMap_noise);

    #prepare windows of particle for scaling
    frequencyMap_mapWindow = calculate_frequency_map(np.zeros((wn_locscale, wn_locscale, wn_locscale)));

    for k in xrange(0, sizeMap[0] - int(wn_locscale), stepSize):
        for j in xrange(0, sizeMap[1] - int(wn_locscale), stepSize):
            for i in xrange(0, sizeMap[2] - int(wn_locscale), stepSize):
                
                emmap_wn = emmap[k: k + wn_locscale, j: j + wn_locscale, i: i + wn_locscale];
                modmap_wn = modmap[k: k + wn_locscale, j: j + wn_locscale, i: i + wn_locscale];
                    
                #do sharpening of the sliding window
                emmap_wn_FFT = np.fft.rfftn(np.copy(emmap_wn));
                modmap_wn_FFT = np.fft.rfftn(np.copy(modmap_wn));

                em_profile, frequencies_map = compute_radial_profile(emmap_wn_FFT, frequencyMap_mapWindow);
                mod_profile, _ = compute_radial_profile(modmap_wn_FFT, frequencyMap_mapWindow);
                scale_factors = compute_scale_factors(em_profile, mod_profile);
                map_b_sharpened, map_b_sharpened_FFT = set_radial_profile(emmap_wn_FFT, scale_factors, frequencies_map, frequencyMap_mapWindow, emmap_wn.shape);
                  
                #do interpolation of sharpening factors
                scale_factors_noise = np.interp(frequencies_noise, frequencies_map, scale_factors);

                #scale noise window with the interpolated scaling factors
                mapNoise_sharpened, mapNoise_sharpened_FFT = set_radial_profile(np.copy(noiseMapFFT), scale_factors_noise, frequencies_noise, frequencyMap_noise, noiseMap.shape);

                #local filtering routines
                if locFilt == True:
                    tmpRes = round(apix/locResMap[k, j, i], 3);
                                                
                    mapNoise_sharpened = lowPassFilter(mapNoise_sharpened_FFT, frequencyMap_noise, tmpRes, noiseMap.shape);
                    map_b_sharpened = lowPassFilter(map_b_sharpened_FFT, frequencyMap_mapWindow, tmpRes, emmap_wn.shape);
                    
                    #calculate noise statistics	  
                    map_noise_sharpened_data = mapNoise_sharpened;	   
                        
                    if ecdfBool:
                        tmpECDF, sampleSort = estimateECDFFromMap(map_noise_sharpened_data, -1, -1);
                        ecdf = np.interp(map_b_sharpened[central_pix, central_pix, central_pix], sampleSort, tmpECDF, left=0.0, right=1.0); 
                    else:
                        ecdf = 0;

                    mean = np.mean(map_noise_sharpened_data);
                    var = np.var(map_noise_sharpened_data);

                    if var < 0.05:
                        var = 0.05;
                        mean = 0.0;
                    if tmpRes == round(apix/100.0, 3):  
                        mean = 0.0;
                        var = 0.0;
                        ecdf = 0;   	 			
                else:    
                    #calculate noise statistics
                    map_noise_sharpened_data = np.copy(mapNoise_sharpened);
                 
                    if ecdfBool:
                        tmpECDF, sampleSort = estimateECDFFromMap(map_noise_sharpened_data, -1, -1);
                        ecdf = np.interp(map_b_sharpened[central_pix, central_pix, central_pix], sampleSort, tmpECDF, left=0.0, right=1.0); 
                    else:
                        ecdf = 0;
                                        
                    mean = np.mean(map_noise_sharpened_data);
                    var = np.var(map_noise_sharpened_data);
                    if var < 0.05:
                        var = 0.05;
                        mean = 0.0;

                #put values back into the the original maps
                halfStep=int((wn_locscale/2.0) - (stepSize/2.0));
                sharpened_map[k + halfStep : k + halfStep + stepSize, j + halfStep : j + halfStep + stepSize, i + halfStep : i + halfStep + stepSize] = np.copy(map_b_sharpened[halfStep:halfStep+stepSize, halfStep:halfStep+stepSize, halfStep:halfStep+stepSize]);
                sharpened_mean_vals[k + halfStep : k + halfStep + stepSize, j + halfStep : j + halfStep + stepSize, i + halfStep : i + halfStep + stepSize] = mean;
                sharpened_var_vals[k + halfStep :  k + halfStep + stepSize, j + halfStep : j + halfStep + stepSize, i + halfStep : i + halfStep + stepSize] = var;
                sharpened_ecdf_vals[k + halfStep : k + halfStep + stepSize, j + halfStep : j + halfStep + stepSize, i + halfStep : i + halfStep + stepSize] = ecdf;
       
    return sharpened_map, sharpened_mean_vals, sharpened_var_vals, sharpened_ecdf_vals;

def get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, wn_locscale, apix, locFilt, locResMap, boxCoord, ecdfBool):

    sharpened_vals = [];
    sharpened_mean_vals = [];
    sharpened_var_vals = [];
    sharpened_ecdf_vals = []; 
       
    central_pix = int(round(wn_locscale / 2.0));
    sizeMap = emmap.shape;
    center = np.array([0.5*sizeMap[0], 0.5*sizeMap[1], 0.5*sizeMap[2]]);
    
    #get the background noise sample
    if boxCoord == 0:
        noiseMap = emmap[int(center[0]-0.5*wn):(int(center[0]-0.5*wn) + wn), int(0.02*wn+wn_locscale):(int(0.02*wn+wn_locscale) + wn), (int(center[2]-0.5*wn)):(int((center[2]-0.5*wn) + wn))];
    else:
        noiseMap = emmap[int(boxCoord[0]-0.5*wn + wn_locscale):(int(boxCoord[0]-0.5*wn + wn_locscale) + wn), int(boxCoord[1]-0.5*wn+ wn_locscale):(int(boxCoord[1]-0.5*wn + wn_locscale) + wn), (int(boxCoord[2]-0.5*wn + wn_locscale)):(int((boxCoord[2]-0.5*wn + wn_locscale)+wn))];

    #prepare noise map for scaling
    frequencyMap_noise = calculate_frequency_map(noiseMap);
    noiseMapFFT = np.fft.rfftn(noiseMap);
    noise_profile, frequencies_noise = compute_radial_profile(noiseMapFFT, frequencyMap_noise);

    #prepare windows of particle for scaling
    frequencyMap_mapWindow = calculate_frequency_map(np.zeros((wn_locscale, wn_locscale, wn_locscale)));

    for k, j, i in (masked_xyz_locs - wn_locscale / 2.0):
		
        emmap_wn = emmap[k: k+wn_locscale, j: j+wn_locscale, i: i+ wn_locscale];
        modmap_wn = modmap[k: k+wn_locscale, j: j+wn_locscale, i: i+ wn_locscale];
	    
        #do sharpening of the sliding window
        emmap_wn_FFT = np.fft.rfftn(np.copy(emmap_wn));
        modmap_wn_FFT = np.fft.rfftn(np.copy(modmap_wn));

        em_profile, frequencies_map = compute_radial_profile(emmap_wn_FFT, frequencyMap_mapWindow);
        mod_profile, _ = compute_radial_profile(modmap_wn_FFT, frequencyMap_mapWindow);
        scale_factors = compute_scale_factors(em_profile, mod_profile);
        map_b_sharpened, map_b_sharpened_FFT = set_radial_profile(emmap_wn_FFT, scale_factors, frequencies_map, frequencyMap_mapWindow);
	  
        #do interpolation of sharpening factors
        scale_factors_noise = np.interp(frequencies_noise, frequencies_map, scale_factors);

        #scale noise window with the interpolated scaling factors
        mapNoise_sharpened, mapNoise_sharpened_FFT = set_radial_profile(np.copy(noiseMapFFT), scale_factors_noise, frequencies_noise, frequencyMap_noise);

        #local filtering routines
        if locFilt == True:
            tmpRes = round(apix/locResMap[k, j, i], 3);
            				
            mapNoise_sharpened = lowPassFilter(mapNoise_sharpened_FFT, frequencyMap_noise, tmpRes);
            map_b_sharpened = lowPassFilter(map_b_sharpened_FFT, frequencyMap_mapWindow, tmpRes);
            
            #calculate noise statistics	  
            map_noise_sharpened_data = mapNoise_sharpened;	   
          	
            if ecdfBool:
                tmpECDF, sampleSort = estimateECDFFromMap(map_noise_sharpened_data, -1, -1);
                ecdf = np.interp(map_b_sharpened[central_pix, central_pix, central_pix], sampleSort, tmpECDF, left=0.0, right=1.0); 
            else:
                ecdf = 0;

            mean = np.mean(map_noise_sharpened_data);
            var = np.var(map_noise_sharpened_data);

            if var < 0.05:
                var = 0.05;
                mean = 0.0;
            if tmpRes == round(apix/100.0, 3):  
                mean = 0.0;
                var = 0.0;
                ecdf = 0;   	 			
        else:    
            #calculate noise statistics
            map_noise_sharpened_data = np.copy(mapNoise_sharpened);
         
            if ecdfBool:
                tmpECDF, sampleSort = estimateECDFFromMap(map_noise_sharpened_data, -1, -1);
                ecdf = np.interp(map_b_sharpened[central_pix, central_pix, central_pix], sampleSort, tmpECDF, left=0.0, right=1.0); 
            else:
                ecdf = 0;
            			
            mean = np.mean(map_noise_sharpened_data);
            var = np.var(map_noise_sharpened_data);
            if var < 0.05:
                var = 0.05;
                mean = 0.0;

        #append values to the sharpened values
        sharpened_vals.append(map_b_sharpened[central_pix, central_pix, central_pix]);
        sharpened_mean_vals.append(mean);
        sharpened_var_vals.append(var);     
        sharpened_ecdf_vals.append(ecdf);
       
    return np.array(sharpened_vals, dtype=np.float32), np.array(sharpened_mean_vals, dtype=np.float32), np.array(sharpened_var_vals, dtype=np.float32), np.array(sharpened_ecdf_vals, dtype=np.float32);

def put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape):
    map_scaled = np.zeros(np.prod(map_shape))
    map_scaled[masked_indices] = sharpened_vals
    map_scaled = map_scaled.reshape(map_shape)
        
    return map_scaled;

def run_window_function_including_scaling(emmap, modmap, mask, wn, wn_locscale, apix, locFilt, locResMap, boxCoord, ecdfBool):
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
    masked_xyz_locs, masked_indices, map_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn_locscale);

    sharpened_vals, sharpened_mean_vals, sharpened_var_vals, sharpened_ecdf_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, wn_locscale, apix, locFilt, locResMap, boxCoord, ecdfBool);
     
    map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape);
    mean_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_mean_vals, masked_indices, map_shape);	 
    var_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_var_vals, masked_indices, map_shape);
    ecdf_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_ecdf_vals, masked_indices, map_shape);

    return map_scaled, mean_map_scaled, var_map_scaled, ecdf_map_scaled; 
   
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
   

def run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, wn_locscale, apix, locFilt, locResMap, boxCoord, ecdfBool):
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
    from mpi4py import MPI;
    comm = MPI.COMM_WORLD;
    rank = comm.Get_rank();
    size = comm.Get_size();
         
    if rank == 0:

        print("*****************************")
        print("********* LocScale **********")
        print("*****************************")

        masked_xyz_locs, masked_indices, map_shape = \
        get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn);
         
        zs, ys, xs = masked_xyz_locs.T;
        zs = split_sequence_evenly(zs, size);
        ys = split_sequence_evenly(ys, size);
        xs = split_sequence_evenly(xs, size);
    else:
        zs = None;
        ys = None;
        xs = None;
     
    zs = comm.scatter(zs, root=0);
    ys = comm.scatter(ys, root=0);
    xs = comm.scatter(xs, root=0);
 
    masked_xyz_locs = np.column_stack((zs, ys, xs));
 
    sharpened_vals, sharpened_mean_vals, sharpened_var_vals, sharpened_ecdf_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, wn_locscale, apix, locFilt, locResMap, boxCoord, ecdfBool);

    comm.barrier();
    sharpened_vals = comm.gather(sharpened_vals, root=0);
    sharpened_mean_vals = comm.gather(sharpened_mean_vals, root=0);
    sharpened_var_vals = comm.gather(sharpened_var_vals, root=0);		
    sharpened_ecdf_vals = comm.gather(sharpened_ecdf_vals, root=0);
     
    if rank == 0:
        sharpened_vals = merge_sequence_of_sequences(sharpened_vals);
        sharpened_mean_vals = merge_sequence_of_sequences(sharpened_mean_vals);
        sharpened_var_vals = merge_sequence_of_sequences(sharpened_var_vals);
        sharpened_ecdf_vals = merge_sequence_of_sequences(sharpened_ecdf_vals);

        map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_vals),
        masked_indices, map_shape);

        mean_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_mean_vals),
        masked_indices, map_shape);

        var_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_var_vals),
        masked_indices, map_shape);
        
        ecdf_map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_ecdf_vals),
        masked_indices, map_shape);

    else:
        map_scaled = None;
        mean_map_scaled = None;
        var_map_scaled = None;		
        ecdf_map_scaled = None;

    comm.barrier();		

    return map_scaled, mean_map_scaled, var_map_scaled, ecdf_map_scaled, rank;
  
def write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol, filename):
    
    if window_bleed_and_pad:
        map_shape = [(LocScaleVol.shape[0] - wn), (LocScaleVol.shape[1] - wn), (LocScaleVol.shape[2] - wn)]
        LocScaleVol = pad_or_crop_volume(LocScaleVol, (map_shape))
                                         
    with mrcfile.new(filename, overwrite=True) as LocScaleVol_out:
        LocScaleVol_out.set_data(LocScaleVol.astype(np.float32))
        LocScaleVol_out.voxel_size = np.rec.array(( args.apix,  args.apix,  args.apix), dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])
        LocScaleVol_out.header.nxstart, LocScaleVol_out.header.nystart, LocScaleVol_out.header.nzstart = [0,0,0]

    return LocScaleVol;

def launch_amplitude_scaling(args):

    startTime = time.time();
    emmap, modmap, mask, wn, wn_locscale, window_bleed_and_pad, method, locFilt, locResMap, boxCoord = prepare_mask_and_maps_for_scaling(args); 
    meanNoise, varNoise, sample = estimateNoiseFromMap(emmap, wn, boxCoord);

    #set output filenames
    if args.outputFilename is not None:
        splitFilename = os.path.splitext(os.path.basename(args.outputFilename))
    else:
        splitFilename = os.path.splitext(os.path.basename(args.em_map))

    if args.testProc is not None:
        testProc = args.testProc
    else:
        testProc = 'rightSided'

    if not args.mpi:
        
        if args.stepSize is None:
            stepSize = 5;
            LocScaleVol, meanVol, varVol, ecdfVol = calculate_scaled_map(emmap, modmap, mask, wn, wn_locscale, args.apix, locFilt, locResMap, boxCoord, args.ecdf, stepSize);
        else:
            stepSize = int(args.stepSize);
            if stepSize == 1:
	            LocScaleVol, meanVol, varVol, ecdfVol = run_window_function_including_scaling(emmap, modmap, mask, wn, wn_locscale , args.apix, locFilt, locResMap, boxCoord, args.ecdf);
            elif stepSize <= 0:
                print("Invalid step size parameter. It has to be greater than 0! Quit program ...");
                return;
            else:
                LocScaleVol, meanVol, varVol, ecdfVol = calculate_scaled_map(emmap, modmap, mask, wn, wn_locscale, args.apix, locFilt, locResMap, boxCoord, args.ecdf, stepSize);
        
        print("Local amplitude scaling finished ...")

        LocScaleVol = mask*LocScaleVol;


        if not args.ecdf:
            ecdfVol = 0;

        qVol = calcQMap(LocScaleVol, meanVol, varVol, ecdfVol, 0, 0, mask, method, testProc);
        qVol = np.subtract(np.ones(qVol.shape), qVol);

        #write the volumes		
        LocScaleVol = write_out_final_volume_window_back_if_required(args, wn_locscale, window_bleed_and_pad, LocScaleVol, splitFilename[0] + '_scaled.mrc');
        qVol = write_out_final_volume_window_back_if_required(args, wn_locscale, window_bleed_and_pad, qVol, splitFilename[0] + '_confidenceMap.mrc');

        endTime = time.time()
        runTime = endTime - startTime
        makeDiagnosticPlot(emmap, wn, wn_locscale, True, boxCoord);
        printSummary(args, runTime)

    elif args.mpi:
        LocScaleVol, meanVol, varVol, ecdfVol, rank = run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, wn_locscale , args.apix, locFilt, locResMap, boxCoord, args.ecdf);       
        if rank == 0:
            print("Local amplitude scaling finished ...")			

            if not args.ecdf:
                ecdfVol = 0;

            qVol = calcQMap(LocScaleVol, meanVol, varVol, ecdfVol, 0, 0, mask, method, testProc);
            qVol = np.subtract(np.ones(qVol.shape), qVol);

            #write the volumes
            LocScaleVol = write_out_final_volume_window_back_if_required(args, wn_locscale, window_bleed_and_pad, LocScaleVol, splitFilename[0] + '_scaled.mrc')
            qVol = write_out_final_volume_window_back_if_required(args, wn_locscale, window_bleed_and_pad, qVol, splitFilename[0] + '_confidenceMap.mrc')

            endTime = time.time()
            runTime = endTime - startTime;
            makeDiagnosticPlot(emmap, wn, wn_locscale, True, boxCoord);
            printSummary(args, runTime);
