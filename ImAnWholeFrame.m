clear all;

%get images
[FileNameMRC, PathNameMRC] = uigetfile('Please select the image file', 'MultiSelect','on');
numMicrographs = length(FileNameMRC);


%create gui and ask for input
prompt = {'FDR:', 'pixel Size', 'min. Size [#Pixel]:', 'binning Factor'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'0.1', '1', '0', '4'};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
FDRval = str2num(answer{1});
pixA = str2num(answer{2});
minSize = str2num(answer{3});
binFac = str2num(answer{4}); 

%ask for directory where output should be saved
outputPath = uigetdir('', 'Please select a folder for the ouput');

%************************************************************************
%*************************** read movies ********************************
%************************************************************************

for micInd = 1:numMicrographs
    %get pure filename
    [ ~, fileName, ~] = fileparts(FileNameMRC{micInd});

    %read first parts of the movie
    [map, s] = ReadMRC(fullfile(PathNameMRC, FileNameMRC{micInd}), 1, 1);
    map = double(map);
    %map = flattenImage(map);

    %define summed image
    mapSum = zeros(s.nx, s.ny);
    mapSum = mapSum + map;

    %define binned map
    mapBinned = zeros(floor(s.nx/binFac), floor(s.ny/binFac), s.nz);
    binnedMapSize = size(mapBinned);

    %normalization
    %mu = mean(mean(map(:, :)));
    %variance = var(reshape(map(:, :), s.nx*s.ny, 1));
    %map( :, :) = (map( :, :) - mu)/sqrt(variance);
    map = localNormalize(map);
    
    
    %do binning
    mapBinned( :, :, 1) = binningImage(map( :, :), binFac);

    %do smoothing with gaussian kernel, i.e. low pass filter
    sigma = 1;
    halfwidth = 3 * sigma;
    [xx, yy] = meshgrid(-halfwidth:halfwidth, -halfwidth:halfwidth);
    smoothMat = exp(-1/(2*sigma^2) * (xx.^2 + yy.^2));
    smoothMat = sum(sum(smoothMat)) * smoothMat;
    mapBinned(:, :, 1) = conv2(mapBinned(:, :, 1), smoothMat, 'same');

    for i = 2:s.nz

        %read frame
        [map, s] = ReadMRC(fullfile(PathNameMRC, FileNameMRC{micInd}), i, 1);
        map = double(map);
        %map = flattenImage(map);

        %add to the average
        mapSum = mapSum + map;	
	
        %normalization
        %mu = mean(mean(map(:, :)));
        %variance = var(reshape(map(:, :), s.nx*s.ny, 1));
        %map( :, :) = (map( :, :) - mu)/sqrt(variance);
        map = localNormalize(map);
        
        %do binning
        mapBinned(:, :, i) = binningImage(map(:, :), binFac );
        clear map; %free memory

        %smoothing
        mapBinned(:, :, i) = conv2(mapBinned(:, :, i), smoothMat, 'same');

    end

    %now perform pixel wise t-tests
    disp('start t-tests');
    pValues = zeros((binnedMapSize(1)*binnedMapSize(2)),1);
    for i = 1:binnedMapSize(1)
        for j = 1:binnedMapSize(2)
            %now calculate t-image, i.e. perform t test for each time series under H0 <= 0
            [~,p] = ttest(mapBinned(i,j,:), 0, 0.05, 'left');	
            pValues((i-1)*binnedMapSize(1)+j) = p;
        end
    end

    disp('start FDR control');
    %do fdr control
    [FDR] = mafdr(pValues, 'BHFDR', true);

    %transform to binary by means of the adjusted p-values
    binImage = zeros(binnedMapSize(1), binnedMapSize(2));
    for i = 1:binnedMapSize(1)
        for j = 1:binnedMapSize(2)
            if( FDR((i-1)*binnedMapSize(1)+j) < FDRval )				
                binImage(i,j) = 1;
            else
                binImage(i,j) = 0;
            end	
        end
    end

    %delete connected components smaller than minsize
    binImage = compAnalysis( binImage, minSize, true);

    mapSum = localNormalize(mapSum); %do normalization for output

    %bring the binary image upon the original size and do point-wise multiplication with the original image
    binImageOriginalSize =  double(imresize(binImage, [s.nx s.ny], 'nearest'));
    outImage = binImageOriginalSize .* mapSum; 


    %now write the mrc files
    WriteMRC(outImage, pixA, fullfile(outputPath, strcat(fileName, '_FDRthresholded.mrc')));  %write thresholded image
    WriteMRC(mapSum, pixA, fullfile(outputPath, strcat(fileName, '_movieSum.mrc')));  %write original movie summed
    WriteMRC(binImageOriginalSize, pixA, fullfile(outputPath, strcat(fileName, '_mask.mrc')));  %write binary mask from fdr thresholding

end
