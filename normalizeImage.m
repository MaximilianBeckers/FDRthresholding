function normalizedIm = normalizeImage( Image ) 
%NORMALIZEIMAGE takes an image as a 2d matrix and normalizes it, where
%the area of the matrix outside a disk with radius radSq is used for
%calculation of noise mean and variance
%
%INPUT: Image: 2d matrix
%Output: normalizedIm: 2d matrix
%
%Author: Maximilian Beckers, EMBL, Structural and Computational Biology, Carsten Sachse Group (2016)
 
	%get image size
	picSize = size(Image);
	
	%squared radius for noise calculation
	radSq = (max(picSize(1)/2, picSize(2)/2) - 5)^2;

	%noise vector
	noiseVec = zeros(ceil(picSize(1)*picSize(2) - pi*radSq), 1);

	%determine center of the image
	centerImCrop = [(1+picSize(1))/2 , (1+picSize(2))/2];

	%get elements in cropped image which are used for calculating noise
	%level (i.e. outside a circle centerend in the cropped image) 	
	pixInd = 1;
	for xInd = 1:picSize(1)
		for yInd = 1:picSize(2)
		        if ((xInd-centerImCrop(1))^2 + (yInd-centerImCrop(2))^2) >= radSq  
		            noiseVec(pixInd) = Image(xInd, yInd);
		            pixInd = pixInd + 1;
			end
		end		     	
	end

	%get mean of noise 
	mu = mean(noiseVec);
		    
	%get variance of noise 
	sigmaSq = var(noiseVec); 

	%now finally do normalization
	normalizedIm = (Image - mu)/ sqrt(sigmaSq);

end
