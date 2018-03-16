function scaledData = coarsegrain(data,tau)
%   scaledData = coarsegrain(data,tau)
%
%   OVERVIEW:  generate the consecutive coarse-grained time series
%   
%   INPUTS:    
%     data  : Vector containing the original seires
%     tau   : scale factor (if tau = 1 returs the original series)
%   OUTPUTS:
%     scaledData : the coarse-grained time series at the scale factor tau
%
%   Written by: Giulia Da Poian <giulia.dap@gmail.com>
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox 
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

if nargin < 2
    error('Wrong number of input arguments');
end

L = length(data);
J = fix(L/tau);

scaledData = zeros(1,J);

if tau < 1
    error('The scale factor must be a positive integer number');
elseif tau == 1
    scaledData = data; % return orignal series
    return
else
    for idx=1:J
        scaledData(idx) = mean(data((idx-1)*tau+1:idx*tau));
    end
end


