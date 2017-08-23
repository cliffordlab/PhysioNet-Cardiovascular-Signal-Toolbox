function xOut = coarsegrain(x)

% output = coarsegrain(inputData, factor, method)
% 
% Overview
%   Coarse-grains a time series by averaging two-item sequences,
%   e.g. new_datum = (i'th odd datum + (i'th + 1) datum) / 2;
%   e.g. x = [1 2 3 5 2 3 5 6 8 9]
%     xOut = [(1+2)/2 (3+5)/2 (2+3)/2 (5+6)/2 (8+9)/2]
%     xOut = [1.5 4.0 2.5 5.5 8.5]
%
% Input    
% 	x - 1D vector of ints or doubles
%
% Output
%   xOut - 1D time series vector after coarse-graining
%
% Dependencies
%
% Reference(s)
%
% Authors
%   Erik Reinertsen <er@gatech.edu>
% 
% Copyright (C) 2017 Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.

n = length(x);

% If vector length is not divisible by 'factor'
% decrement last index until it is
lastIndex = n;
while(mod(lastIndex, 2) > 0)
    lastIndex = lastIndex - 1;
end

% Trim inputData until vector length is divisible by 0
x = x(1:lastIndex);

% Calculate length of new inputData vector
lenInputDataNew = floor(n * 0.5 ^ (2 - 1));

% % Initialize vector to store coarse-grained vector
% output = NaN(lenInputDataNew, 1);

% Loop through each index of output
% compute mean of the i'th to (i+factor)th indices,
% and write to i'th element of output
for idxNew = 1:lenInputDataNew
    index = 1 + 2 * (idxNew - 1);
    xOut(idxNew) = mean(x(index:index + 2 - 1));
end

end