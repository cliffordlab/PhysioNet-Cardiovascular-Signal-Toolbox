function output = runrmssd(rr)
% output = runrmssd(times, rr);
%   OVERVIEW:   calculates the rmssd: the square root of the mean 
%               of the sum of the squares of the differences.
%   INPUT:      rr = a single row of rr interval data in seconds
%               times =  the time indices of the rr interval data (seconds)
%   OUTPUT:     output = rmssd
%   DEPENDENCIES
%   & LIBRARIES:    replaceOutliers.m
%                   AnnotationConversion.m
%   REFERENCE:  
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   SOURCE:     Adapted from Kaplan Tools
%               See: http://www.robots.ox.ac.uk/~gari/CODE/HRV/
%   AUTHORS:    Adriana Nicholson Vest, PhD
%   COPYRIGHT (C) 2016
%   LICENSE:    
%%

leng = length(rr);

dif = diff(rr);

output = sqrt(mean(dif.*dif));

end