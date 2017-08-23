function mse = multiscaleEntropy(data, m, r, maxTau, entropyType)

% mse = multiscaleEntropy(data, m, r, max_tau, entropy_script)
%
% Overview
%	Calculates multiscale entropy on a vector of input data.
%
% Input
%	data        - data to analyze; vector of doubles
%   m           - pattern length; int
%   r           - radius of similarity (% of std); double
%   maxTau      - maximum number of coarse-grainings; int
%   entropyType - type of entropy to calculate
%
% Output
%   mse          - vector of [max_tau, 1] doubles
%
% Example
%   data = rand(1e4, 1); % generate random data
%   m = 2; % template length
%   r = 0.2; % radius of similarity
%   entropyType = 'sampen'; % sample entropy, but can also choose 'fuzzy', 
%   and this can be modularly extended for other entropy types
%   maxTau = 4; % calculate sample entropy over four coarse grainings
%   mse = multiscaleEntropy(data, m, r, maxTau, entropyType);
% 
% Dependencies
%   https://github.com/cliffordlab/
%
% Reference(s)
% 
% Copyright (C) 2017 Erik Reinertsen <er@gatech.edu>
% All rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.


% Initialize output vector
mse = NaN(maxTau, 1);

% Loop through each timescale
% Note: i_tau == 1 is the original time series
for i_tau = 1:maxTau
    
    % Calculate entropy
    if strcmp(entropyType, 'sampen')
        mse(i_tau) = fastSampen(data, m, r);
    elseif strcmp(entropyType, 'fuzzy')
        [mse(i_tau), ~] = fuzzyEntropy(data, m, r);
    elseif strcmp(entropyType, 'shannon')
        mse(i_tau) = shannonEntropy(data);
    elseif strcmp(entropyType, 'debug')
        mse(i_tau) = NaN;
    end
    
    % Coarse-grain data (default: halving method)
    data = coarsegrain(data);
    
end % end for loop

end % end function
