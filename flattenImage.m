function flatImage = flattenImage(Image)

	
	%order of polynomial to fit
	n = 2;

	%get image size
	imSize = size(Image);
	flatImage = zeros(imSize(1), imSize(2));
    
	%do flattening for rows
	for xInd = 1:imSize(1)
		xLine = Image(xInd, :);
		indices = (1:imSize(2));
		polyCoeff = polyfit(indices, xLine, n);
		flatImage(xInd, :) = xLine - polyval(polyCoeff, indices);
	end

	%do flattening for columns
	for yInd = 1:imSize(2)
		yLine = flatImage(:, yInd)';
		indices = (1:imSize(1));
		polyCoeff = polyfit(indices, yLine, n);
		flatImage(:, yInd) = yLine - polyval(polyCoeff, indices);
	end	


end
