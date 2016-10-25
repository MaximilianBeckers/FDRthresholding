function Image = removeOutliers( Image )
%REMOVEOUTLIERS takes an image as a matrix and sets the pixel values of outliers to
%the mean of the image. Outliers are detected as being 3*stand.dev. higher or lower
%than the mean of the image.
%
%Input: Image: 2d matrix
% 	
%Output: Image: 2d matrix with removed Outliers
%
%Author: Maximilian Beckers, EMBL, Structural and Computational Biology, Carsten Sachse Group (2016)

	%get image size
	picSize = size(Image);
	
	%get mean, variance and standard deviation from image
	meanImage = mean(mean(Image));
	varImage = var(reshape(Image, picSize(1)*picSize(2), 1));
	sdImage = sqrt(varImage);

	%set upper and lower bounds for pixel values
	upBound = meanImage + 3*sdImage;
	lowBound = meanImage - 3*sdImage;

	%find the outliers according to the bounds and set their values to the mean
	outliers = find((Image > upBound) && (Image < lowBound))
	if ~isempty(outliers)	
		Image(outliers) = meanImage;
	end

end
