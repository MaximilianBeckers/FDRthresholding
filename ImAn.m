%clear workspace
clear all;


%get images
[FileNameMRC,PathNameMRC] = uigetfile('Please select the image files','MultiSelect','on');

%get the boxcoordinates
[FileName,PathName] = uigetfile('Please select the box coordinate files','MultiSelect','on');

%load the box coordinates of all boxes in the list
fileID = fopen(fullfile(PathName,FileName),'r');
boxes = fgetl(fileID, %f	%f	%f	%f);
fclose(fileID);

%transpose to macht original orientation of the .box - file
boxes = boxes'; % x y width height


%**************************************************************************
%********************* read images ****************************************
%**************************************************************************


%read first parts of the image
[map,s]=ReadMRC(fullfile(PathNameMRC{mrcInd}, FileNameMRC{mrcInd}), 1,15);

%binning factor
binFac = 4;

%init new image for binning
mapNew = zeros( round(s.ny/binFac), round(s.nx/binFac) );

%bin and average the image for proper representation
i = 1;
j = 1;
k = 1;
l = 1;

while i < s.nx
    while j < s.ny
        %binning and summing over all frames
        mapNew(l,k) = sum(sum(sum(map(i:(i+(binFac-1)), j:(j+(binFac-1)),1:s.nz))))/(binFac*binFac*s.nz); 
        
        j = j + binFac;
        l = l + 1;
    end  
    j = 1;
    l = 1;
    
    i = i + binFac;
    k = k + 1;
end

%let user select subimage
fig = figure;
imshow(mapNew, 'DisplayRange', []);
rect = getrect(fig);
close(fig);

%transform rect to full image coordinates
col = (rect(1)-1)*binFac + 1;
row = (rect(2)-1)*binFac + 1;
width = round(rect(3)*binFac);
height = round(rect(4)*binFac);
pictSize = [height, width];

%crop in full image
cropDat = map(col:(col+width-1), row:(row+height-1), :);

%radius for noise calculation
rad = (max(width/2, height/2) - 1)^2;

%normalize and flatten images
for i = 1:s.nz
    
    centerImCrop = [(1+height)/2 , (1+width)/2];

    %flatten image
    %first lines
    %for j = 1: pictSize(1)
    %   indexVec = (1:pictSize(2));  
    %   vec = cropDat(j,:,i);
    %   p = polyfit(indexVec,vec,n);
    %   polVal = polyval(p,indexVec);
    %   cropDat(j, :, i) = cropDat(j, :, i) - polVal;
    %end

    %then columns
    %for j = 1: pictSize(2)
    %   indexVec = (1:pictSize(1));  
    %   vec = cropDat(:,j, i)';
    %   p = polyfit(indexVec,vec,n);
    %   polVal = polyval(p,indexVec);
    %   cropDat(:, j, i) = cropDat(:, j, i)' - polVal;
    %end
    
    %get elements in cropped image which are used for calculating noise level (i.e. outside a circle centerend in the cropped image) 	
    noiseVec = [];
    for rowInd = 1:pictSize(1)
        for colInd = 1:pictSize(2)
            if (rowInd-centerImCrop(1))^2 + (colInd-centerImCrop(2))^2 >= rad  
                noiseVec = [noiseVec, cropDat(colInd, rowInd, i)];	
            end		     	
        end
    end 

    %get mean of noise in frame i
    mu = mean(noiseVec);
    
    %get variance of noise in frame i
    sigmaSq = var(noiseVec); 

    %now finally do normafglization
    cropDat(:,:, i) = (cropDat(:,:, i) - mu)/ sqrt(sigmaSq);
    
    %do ctf correction
    %fm = fftshift(fft2(cropDat(:, :, i))); %fourier transformation of image
    %parameters
    pixA = 1.77;
    lambda = 0.0197;
    Cs = 2.7;
    defocus = 9.28632;
    B = 0;
    alpha = 0.07;
    theta = 25.50127;
    deltadef = 0.05668;
    %multipl. in fourier space and inverse fourier trafo
    %cropDat(: , :, i) = ifft2(ifftshift(fm.*sign(ctf(length(cropDat(:, :, i)), pixA, lambda, defocus, Cs, B, alpha, deltadef, theta)))); 
    
    %do smoothing with gaussian kernel, i.e. low pass filter
    sigma = 2;
    halfwidth = 3 * sigma;
    [xx, yy] = meshgrid(-halfwidth:halfwidth, -halfwidth:halfwidth);
    smoothMat = exp(-1/(2*sigma^2) * (xx.^2 + yy.^2));
    smoothMat = sum(sum(smoothMat)) * smoothMat;
    cropDat(:, :, i) = conv2(cropDat(:, :, i), smoothMat, 'same');

end

clear mapNew;

%now calculate t-image, i.e. perform t test for each time series under H0 = 0
pValues = zeros((height*width),1);
for i = 1:height
	for j = 1:width
        %[~,p] = ttest(cropDat(j,i,:));

		[~,p] = ttest(cropDat(j,i,:),0,0.05, 'left');	
        pValues((i-1)*width+j) = p;
	end
end

%do benjamini-hochberg fdr control
[FDR, q] = mafdr(pValues);

%transform to binary by means of the adjusted p-values
binImage = zeros(height, width);
for i = 1:height
	for j = 1:width
		if( q((i-1)*width+j) < 0.20 )
			%cropIm(i,j) = pValues((i-1)*width+j);
            binImage(i,j) = 1;
		else
			binImage(i,j) = 0;
		end	
	end
end

%delete components that are so small to be neglected
CC = bwconncomp(binImage); % find conn. components
numObj = CC.NumObjects; %number of all components
numPixels = CC.PixelIdxList; %number of pixels for each component
minSize = 20; %minimal number of pixels in component

%now erase all components smaller than minsize
for compInd = 1:numObj
    if length(numPixels{compInd}) < minSize
        binImage(CC.PixelIdxList{compInd}) = 0; %erase
    end
end    

origIm = (imcrop(sum(map,3), [row, col, height, width]))';

figure;
subplot(2,2,1);
imshow(binImage);
subplot(2,2,2);
imshow(origIm, 'DisplayRange', []);
subplot(2,2,3);
imshow( conv2(origIm, smoothMat, 'same'), 'DisplayRange', [] );
