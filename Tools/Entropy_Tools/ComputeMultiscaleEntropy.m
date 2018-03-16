 function mse = ComputeMultiscaleEntropy(data, m, r, maxTau, maxVecSize)

% mse = ComputeMultiscaleEntropy(data, m, r, max_tau)
%
% Overview
%	Calculates multiscale entropy on a vector of input data.
%
% Input
%	data        - data to analyze; vector of doubles
%   m           - pattern length; int
%   r           - radius of similarity (% of std); double
%   maxTau      - maximum number of coarse-grainings; int
%   maxVecSize  - optional, parameter to switch from SampEn to FastSempEn        
%                
% Output
%   mse          - vector of [max_tau, 1] doubles
%
% Example
%   data = rand(1e4, 1); % generate random data
%   m = 2; % template length
%   r = 0.2; % radius of similarity
%   maxTau = 4; % calculate sample entropy over four coarse grainings
%   mse = multiscaleEntropy(data, m, r, maxTau, entropyType);
% 
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
% Reference(s)
% 
% Copyright (C) 2017 Erik Reinertsen <er@gatech.edu>
% All rights reserved.
%
% 
% 08-23-2017 Modyfied by Giulia Da Poian for the VOSIM HRV Toolbox 
% Removed the possibility to use different types of entropy, only
% fastSampen method in this version
%
% 10-10-2017 Modyfied by Giulia Da Poian, use FastSampEn in series < 34000
% otherwise use traditional SampEn that is faster for long series
%
% 10-23-2017 Modyfied by Giulia Da Poian, using scales like in Costa's
% paper instead of Coarse-grain data using halving method.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.

if nargin<5 || isempty(maxVecSize)
    maxVecSize = 34000;
end

data = zscore(data);  % (introduced by GDP) normalization of the signal that 
                      % replace the common practice of expressing the 
                      % tolerance as r times the standard deviation
                      
                      
mse = NaN(maxTau, 1); % Initialize output vector

SampEnType = 'Maxim'; % Initialize default SampEn method 

% Check data length, if < 34000 use Fast Implementation (introduced GDP) 
% if length(data) < maxVecSize
%     SampEnType = 'Fast'; 
% end


% Loop through each window

% Loop through each timescale
% Note: i_tau == 1 is the original time series
for i_tau = 1:maxTau
    % Coarse-grain data (using scale in Costa's paper)
    scaledData = coarsegrain(data,i_tau); % Changed by GDP
    
    switch SampEnType
        case 'Fast'
           mse(i_tau) = fastSampen(scaledData, m, r);
        otherwise
           mse(i_tau) = sampenMaxim(scaledData, m, r); 
    end

end % end for loop

end % end function
