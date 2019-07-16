function ApEn = ApproxEntropy( data, m, r)
% OVERVIEW: This function calculates the approximate entropy for the input
% vector data
%
% INPUTS
%   data : time-series data
%   m    : Template length
%   r    : tolerance (typically 0.2 * std)
% OUTPUT
%   ApEn : approximate entropy value
%
% Written by Giulia Da Poian <giulia.dap@gmail.com>
%       Last modified by Camilo Valderrama -- 7/16/2019
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

N = length(data);
result = zeros(1,2);

for j = 1:2
    m = m+j-1;
    phi = zeros(1,N-m+1);
    dataMat = zeros(m,N-m+1);
    
    % setting up data matrix
    for i = 1:m
        dataMat(i,:) = data(i:N-m+i);
    end
    
    % counting similar patterns using distance calculation
    for i = 1:N-m+1
        
        if sum( isnan( dataMat(:,i) ) ) == 0
        
            tempMat = abs(dataMat - repmat(dataMat(:,i),1,N-m+1));
            boolMat = any( (tempMat > r),1);
            phi(i) = sum(~boolMat)/(N-m+1);
        else
	    %discarding blocks with nan values
            phi(i) = nan;
    
        end
    end
    
    % summing over the counts
    result(j) = log( nanmean( phi) );    
end

ApEn = result(1) - result(2) ;

end