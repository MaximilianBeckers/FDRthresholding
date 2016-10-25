function testResult =  BenjYekut(p, alpha)
%perform Benjamini-Yekutieli procedure

	%number of hypotheses
	m = length(p);

	%output
	testResult = zeros(1,m);

	%sort the p-values
	[pSort, index] = sort(p, 'descend');
	
	%calculate Benj-Yek constant
	c = 0;
	for i = 1:m
		c = c + 1/i; 
	end

	%for simple benjamini-hochberg
	c = 1; 
	
	%perform the step-up procedure
	for i = 1:m 
		if(pSort(i) <= ((m-i+1)/(m*c))*alpha) %reject null hypotheses
			pSort(i:end) = ones(1, length(pSort(i:end)));
			break;						
		else
			pSort(i) = 0;
		end	
	end

	%resort the test results
	for i = 1:m
		testResult(index(i)) = pSort(i); 
	end

end 
