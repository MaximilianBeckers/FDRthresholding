function [binaryIm] = compAnalysis( binaryIm, minSize, wholeFrame)
	
	%get image size
	pictSize = size(binaryIm);
	
	maxRadsq = (max(pictSize(1), pictSize(2))*0.5)^2;
	
	%get center of image
	center = [pictSize(1)/2, pictSize(2)/2];
	
	%do connected component analysis		
	CC = bwconncomp(binaryIm);
	numComp = CC.NumObjects;
	
	%now check each component if it is outside the ball	
	for compInd = 1:numComp
	
		%get pixels in the respective component
		pixelList = CC.PixelIdxList{compInd};
        
        if length(pixelList) < minSize
            binaryIm(CC.PixelIdxList{compInd}) = 0; %erase
            continue;
	    end
		
		if ~wholeFrame
			%transform linear indices to usual matrix indices
			[I,J] = ind2sub(pictSize, pixelList); 
		
			%calculate distances for all pixels to the image center
			dist = (I-center(1)).^2 + (J-center(2)).^2; 
			
			%finally we check if components are in the forbidden region
			if isempty(find(dist > maxRadsq))
				%do nothing, all pixels of this component are in the allowed range
			else
				%delete all pixels of the respective component
				binaryIm(CC.PixelIdxList{compInd}) = 0;
			end
		end
	end
end


