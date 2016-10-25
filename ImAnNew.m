%clear workspace
clear all;

%get movie frames
[tempFileNameMRC, tempPathNameMRC] = uigetfile('*.mrcs','Please select the image files','MultiSelect', 'on');

%get box coordinates
[tempFileName, tempPathName] = uigetfile('*.box','Please select the box coordinate files', 'MultiSelect', 'on');

%get the CTF files
[tempFileNameCTF, tempPathNameCTF] = uigetfile('Please select the ctf data file', 'MultiSelect', 'on');

numMicrographs = length(tempFileNameMRC);
numParticlesPerMic = zeros(numMicrographs, 1);

FileNameMRC  = cell(numMicrographs, 1);
FileName = cell(numMicrographs, 1);
FileNameCTF = cell(numMicrographs, 1);
outNameFDR = cell(numMicrographs, 1);
outNameRaw = cell(numMicrographs, 1); 

%now look for matching image and box and ctf files and do the matching
for micInd = 1:numMicrographs

	FileNameMRC{micInd} = tempFileNameMRC{micInd};
	PathNameMRC= tempPathNameMRC;

	%search the boxcoordinates for each micrograph, this is necessary as the files do not have to be in the same order

	for micInd2 = 1:numMicrographs
		FileName{micInd} = tempFileName{micInd2};
		PathName = tempPathName;
		
		%get the names without paths and file extensions
		[~, fileNameBox, ~] = fileparts(FileName{micInd});
		[~, fileNameMRC, ~] = fileparts(FileNameMRC{micInd});
		
		%check if the filenames are the same at the beginning of the string
		if strncmpi(fileNameBox, fileNameMRC, min(length(fileNameBox), length(fileNameMRC)))
			
			%do the same search now for ctf files
			for micInd3 = 1:numMicrographs

				FileNameCTF{micInd} = tempFileNameCTF{micInd3};
				PathNameCTF= tempPathNameCTF;
				
				%get the names without paths and file extensions
				[~, fileNameCTF, ~] = fileparts(FileNameCTF{micInd});
				
				%again check if the filenames are the same at the beginning
				if strncmpi(fileNameCTF, fileNameMRC, min(length(fileNameCTF), length(fileNameMRC)))	
					break; %obviously we have a match
				end
					
			end
			break; %obviously we have a match	
		end
	end	
end



%*************************************************************************
%************************** gather user input ****************************
%*************************************************************************

%create gui and ask for input
prompt = {'Enter FDR:', 'PixelSize [A]:', 'min. Size [#Pixel]:', 'binFac'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'0.1', '1', '10', '2'};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
FDRval = str2num(answer{1});
pixA = str2num(answer{2});
minSize = str2num(answer{3});
binFac = str2num(answer{4});

%ask for directory where output should be saved
outputPath = uigetdir('','Please select a folder for the ouput');


%*************************************************************************
%************************ analyse all the mics ***************************
%*************************************************************************

for micInd = 1:numMicrographs

	%load the box coordinates of all boxes in the list
	fileID = fopen(fullfile(PathName, FileName{micInd}),'r');
	boxes = fscanf(fileID, '%f %f %f %f', [4,inf] );
	fclose(fileID);

	%transpose to macht original orientation of the .box - file
	boxes = boxes'; % x y width height
	width = boxes(1, 3);
	height = boxes(1, 4);

	[numParticles, ~] = size(boxes);
	numParticlesPerMic(micInd) = numParticles;

	particles = cell(numParticles,1); %define particle cell
	particlesOut = zeros(width, height, numParticles); %initialize output of signif-masked particles
	particlesOutRaw = zeros(width, height, numParticles); %initialize output of raw particles

	%************************************************************************
	%*************************** read movies ********************************
	%************************************************************************

	%read first frame of the movie
	[map, s] = ReadMRC(fullfile(PathNameMRC, FileNameMRC{micInd}), 1,1);

	%now append all particles
	for particleInd = 1:numParticles
		xmin = boxes(particleInd,1);
	    width = boxes(particleInd,3);
		height = boxes(particleInd,4);
	    ymin = boxes(particleInd,2);
	
		%crop the particles in the original map and invert contrast
		particles{particleInd} = imcomplement(map(xmin:(xmin+width-1), ymin:(ymin+height-1)));	
	end

	%free memory for new image
	clear map;

	%now read all the other frames
	for frameInd = 2:s.nz
		[map,s]=ReadMRC(fullfile(PathNameMRC, FileNameMRC{micInd}), frameInd, 1);

		%now append the particles
		for particleInd = 1:numParticles
			xmin = boxes(particleInd,1);
		    width = boxes(particleInd,3);
			height = boxes(particleInd,4);
			ymin = boxes(particleInd,2);
			%crop the particles in the original map

			particles{particleInd} = cat(3, particles{particleInd}, imcomplement(map(xmin:(xmin+width-1), ymin:(ymin+height-1))));	
		end

		%free memory
		clear map;
	end


	%**************************************************************************
	%******************** now do analysis for each particle *******************
	%**************************************************************************

	binnedSize = [floor((width/binFac)), floor((height/binFac))];
	binnedDat = zeros(binnedSize(1), binnedSize(2), s.nz);

	%radius for noise calculation
	rad = (max(width/2, height/2) - 5)^2;

	%noise vector
	noiseVec = zeros(ceil(height*width - pi*rad), 1);
	for particleInd = 1:numParticles

		%crop the right particle
		cropDat = particles{particleInd};

		%normalize frames
		for i = 1:s.nz
		    
		    centerImCrop = [(1+height)/2 , (1+width)/2];
		    %get elements in cropped image which are used for calculating noise
		    %level (i.e. outside a circle centerend in the cropped image) 	
            pixInd = 1;
            for rowInd = 1:height
                for colInd = 1:width
                    if ((rowInd-centerImCrop(1))^2 + (colInd-centerImCrop(2))^2) >= rad  
                        noiseVec(pixInd) = cropDat(colInd, rowInd, i);
                        pixInd = pixInd + 1;
                    end		     	
                end
            end 

		    %get mean of noise in frame i
		    mu = mean(noiseVec);
		    
		    %get variance of noise in frame i
		    sigmaSq = var(noiseVec); 

		    %now finally do normalization
		    cropDat(:,:, i) = (cropDat(:,:, i) - mu)/ sqrt(sigmaSq);
		
            %do binning
            binnedDat(:,:, i) = binningImage(squeeze(cropDat(:,:,i)), binFac);

            %do smoothing with gaussian kernel, i.e. low pass filter
		    sigma = 1;
		    halfwidth = 3 * sigma;
		    [xx, yy] = meshgrid(-halfwidth:halfwidth, -halfwidth:halfwidth);
		    smoothMat = exp(-1/(2*sigma^2) * (xx.^2 + yy.^2));
		    smoothMat = sum(sum(smoothMat)) * smoothMat;
		    binnedDat(:, :, i) = conv2(binnedDat(:, :, i), smoothMat, 'same');
		 
	    end
	    
	    %push back the normalized particles
	    particles{particleInd} = cropDat;

	    disp('start t-tests');
		%now calculate t-image, i.e. perform t test for each time series under H0 = 0
		pValues = zeros((binnedSize(1)*binnedSize(2)),1);
		for i = 1:binnedSize(2)
			for j = 1:binnedSize(1)
				[~,p] = ttest(binnedDat(j,i,:), 0, 0.05, 'right');	
				pValues((i-1)*width+j) = p;
			end
	    end
	    
	    disp('start FDR control');
		%do fdr control
		[FDR] = mafdr(pValues, 'BHFDR', true);

	    disp('start binary conversion');
		%transform to binary by means of the adjusted p-values
		binImage = zeros(binnedSize(2), binnedSize(1));
		for i = 1:binnedSize(2)
			for j = 1:binnedSize(1)
				if( FDR((i-1)*width+j) < FDRval )				
                    			binImage(i,j) = 1;
				else
					binImage(i,j) = 0;
				end	
			end
        end

        %delete small components and components obviously not a particle
	    disp('delete small components');
	    binImage = compAnalysis( binImage, minSize, false );
	    
		%get summed original image for visualization and normalize it
		origIm = (sum(particles{particleInd},3))';
		particles{particleInd} = []; % just free memory
		origIm = normalizeImage(origIm);

		%get output image, i.e. the masked original particle and append it to output
		binImage = imresize(binImage, [width height] );
		outIm = origIm.*binImage;
		particlesOut(:, :, particleInd) =  outIm';
		particlesOutRaw(:, :, particleInd) =  origIm';
	
		%do plotting
		%figure;
		%subplot(2,2,1);
		%imshow(binImage);
		%subplot(2,2,2);
		%imshow(origIm, 'DisplayRange', []);
		%subplot(2,2,3);
		%imshow( conv2(binningImage(origIm, binFac), smoothMat, 'same'), 'DisplayRange', [] );
		%subplot(2,2,4);
		%imshow(outIm, 'DisplayRange', []);
	end


	%************************************************************************
	%*************************** create outputs *****************************
	%************************************************************************
	
	outNameFDR{micInd} = strcat('imstack_', FileNameMRC{micInd});
	outNameRaw{micInd} = strcat('imstackRAW_', FileNameMRC{micInd});

	%create full output filename
	fullOutNameFDR = fullfile(outputPath, outNameFDR{micInd});
	fullOutNameRaw = fullfile(outputPath, outNameRaw{micInd});

	%write MRC output stack
	WriteMRC(particlesOut(:, :, :), pixA, fullOutNameFDR); 
	WriteMRC(particlesOutRaw(:, :, :), pixA, fullOutNameRaw); 

end


%create a new Directory for the images
%newfolder = fullfile(outputPath, 'images'); 
%mkdir(newfolder);

%get the total number of particles
%[~, ~, numAllParticles] = size(particlesOut);

%finally write the images
%for particleInd = 2:numAllParticles
%	imageName = strcat('pic', num2str(particleInd), '.png'); %to do: give names accrodig to micrograph
%	imOutName = fullfile(newfolder, imageName );
%	imwrite(squeeze(particlesOut(:,:, particleInd)), imOutName);
%end


%*************************************************************************
%********************* print .star file for relion ***********************
%*************************************************************************

%print significance masked data

%print header
fileID = fopen(fullfile( outputPath,'signifMaskedImages.star'),'w'); 
fprintf(fileID,'data_images \n');
fprintf(fileID,'loop_ \n');
fprintf(fileID,'_rlnImageName \n');
fprintf(fileID,'_rlnMicrographName \n');
fprintf(fileID,'_rlnDefocusU \n');
fprintf(fileID,'_rlnDefocusV \n');
fprintf(fileID,'_rlnDefocusAngle \n');
fprintf(fileID,'_rlnVoltage \n');
fprintf(fileID,'_rlnSphericalAberration \n');
fprintf(fileID,'_rlnAmplitudeContrast \n');

%print data

for micInd = 1:numMicrographs
	tempFileID = fopen(fullfile(PathNameCTF, FileNameCTF{micInd}),'r');
	ctfDat = fscanf(tempFileID, '%f %f %f ', [3,inf] );
	fclose(tempFileID);
	ctfDat = ctfDat';

    for particleInMic = 1:numParticlesPerMic(micInd)
        fprintf(fileID, '%u@%s %s %f %f %f %f %f %f \n', particleInMic, outNameFDR{micInd}, outNameFDR{micInd}, ctfDat(1, 1), ctfDat(1, 2), ctfDat(1, 3), 300, 2.7, 0.07);
    end
end


%print raw data

%print header
fileID = fopen(fullfile( outputPath, 'images.star'),'w'); 
fprintf(fileID,'data_images \n');
fprintf(fileID,'loop_ \n');
fprintf(fileID,'_rlnImageName \n');
fprintf(fileID,'_rlnMicrographName \n');
fprintf(fileID,'_rlnDefocusU \n');
fprintf(fileID,'_rlnDefocusV \n');
fprintf(fileID,'_rlnDefocusAngle \n');
fprintf(fileID,'_rlnVoltage \n');
fprintf(fileID,'_rlnSphericalAberration \n');
fprintf(fileID,'_rlnAmplitudeContrast \n');

%print data

for micInd = 1:numMicrographs
	tempFileID = fopen(fullfile(PathNameCTF, FileNameCTF{micInd}),'r');
	ctfDat = fscanf(tempFileID, '%f %f %f ', [3,inf] );
	fclose(tempFileID);
	ctfDat = ctfDat';

    for particleInMic = 1:numParticlesPerMic(micInd)
        fprintf(fileID, '%u@%s %s %f %f %f %f %f %f \n', particleInMic, outNameRaw{micInd}, outNameRaw{micInd}, ctfDat(1, 1), ctfDat(1, 2), ctfDat(1, 3), 300, 2.7, 0.07);
    end
end

fclose(fileID);

