function Image = localNormalize( Image )
    
	pictSize = size(Image);
	subImSize = [1024, 1024];

	numXSegment = ceil(pictSize(1)/subImSize(1));
	numYSegment = ceil(pictSize(2)/subImSize(2));

	for xSegInd = 1:numXSegment
		for ySegInd = 1:numYSegment

			%get the respective start and end indices of subimages to normalize at whole
			if xSegInd ~= numXSegment && ySegInd ~= numYSegment

				xStart = (xSegInd-1)*subImSize(1) + 1;
				xEnd = (xSegInd*subImSize(1));
				yStart = ((ySegInd-1)*subImSize(2) + 1);
				yEnd = (ySegInd*subImSize(2));

			elseif xSegInd == numXSegment && ySegInd ~= numYSegment
				
				xStart = (xSegInd-1)*subImSize(1) + 1;
				xEnd = pictSize(1);
				yStart = ((ySegInd-1)*subImSize(2) + 1);
				yEnd = (ySegInd*subImSize(2));
					
 			elseif xSegInd ~= numXSegment && ySegInd == numYSegment
				
				xStart = (xSegInd-1)*subImSize(1) + 1;
				xEnd = (xSegInd*subImSize(1));
				yStart = ((ySegInd-1)*subImSize(2) + 1);
				yEnd = pictSize(2);
			else

				xStart = (xSegInd-1)*subImSize(1) + 1;
				xEnd = pictSize(1);
				yStart = ((ySegInd-1)*subImSize(2) + 1);
				yEnd = pictSize(2);
			end		
 		
			%get subimage and calc mean and var
			subImage = Image( xStart:xEnd, yStart:yEnd);
			meanSubIm = mean(mean(subImage));
			sizeSubIm = size(subImage);						
			varSubIm = var(reshape(subImage, sizeSubIm(1)*sizeSubIm(2), 1));
			
			%do the actual normalization
			subImage = (subImage - meanSubIm)/ sqrt(varSubIm);
			Image( xStart:xEnd, yStart:yEnd) = subImage;		
	
		end
	end
end
