function [binnedIm] = binningImage(Image, binFac)
% This function takes an image and bins it according to the specified binning factor
%
% Input: Image: Image to bin
%	 binFac: number of pixels in x and y direction to sum up
%
% Output: binnedIm: binned Image
%
% Author: Maximilian Beckers, 2016, EMBL

	pictSize = size(Image);
	
	binnedIm = zeros(floor(pictSize(1)/binFac), floor(pictSize(2)/binFac));	
	binnedPictSize = [floor(pictSize(1)/binFac), floor(pictSize(2)/binFac)];
	
	%do the binning
	k = 1;
	l = 1;

	for i = 1:binnedPictSize(1)
		for j = 1:binnedPictSize(2)
			if j == binnedPictSize(2) && i ~= binnedPictSize(1)
				binnedIm(i,j) = sum(sum(Image(k:(k+binFac-1), l:end)));
			elseif i == binnedPictSize(1) && j ~= binnedPictSize(2)
				binnedIm(i,j) = sum(sum(Image(k:end, l:(l+binFac-1))));
			elseif j == binnedPictSize(2) && i == binnedPictSize(1)
				binnedIm(i,j) = sum(sum(Image(k:end, l:end)));
			else
				binnedIm(i,j) = sum(sum(Image(k:(k+binFac-1), l:(l+binFac-1))));
			end			
			l = l + binFac; 	
        end
        l = 1;
		k = k + binFac;
	end
end 
