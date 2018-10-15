 function mse = ComputeMultiscaleEntropy(data, m, r, maxTau, maxVecSize,cg_moment,cg_method,constant_r)

% mse = ComputeMultiscaleEntropy(data, m, r, max_tau,cg_moment)
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
%   cg_moment   - optional, moment used to coarse-grain the time series,
%                 by default use the mean 'mean'
%                 Options: 'mean', 'varaince'
%   cg_method   - [optional] - method use to generate the coarse-grain 
%                 time series. By default use the original method proposed
%                 by Costa et al. ('fir' with mean).
%                 Options: 'fir' [Costa et al.], 'butter' [Valencia et al.]
%   constant_r  - [DEFAULT] 1 : use r as function of std of original time 
%                 series in this implementation correspond to zscore only 
%                 the original time series 
%               - 0 : r as a function of scale factor (tau), i.e., 
%                 r(i)= r*std(ScaleData(i)), in this implementation it
%                 corresponds to zscore  each coarse-grain time siries
% 
%      
% Output
%   mse          - vector of [max_tau, 1] doubles
%
% Example
%   data = rand(1e4, 1); % generate random data
%   m = 2; % template length
%   r = 0.2; % radius of similarity
%   maxTau = 4; % calculate sample entropy over four coarse grainings
%   mse = multiscaleEntropy(data, m, r, maxTau, 'Fast','mean');
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
% 08-23-2017 Modyfied by Giulia Da Poian to be inclused in the HRV Toolbox 
% for PhysioNet Cardiovascular Signal ToolboxToolbox 
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
if nargin<6 || isempty(cg_moment)
    cg_moment = 'mean';
end
if nargin<7 || isempty(cg_method)
    cg_method = 'fir';
end
if nargin<8 || isempty(constant_r)
    constant_r = 0;
end


data = zscore(data);  % (introduced by GDP) normalization of the signal that 
                      % replace the common practice of expressing the 
                      % tolerance as r times the standard deviation
                      
mse = NaN(maxTau, 1); % Initialize output vector

SampEnType = 'Maxim'; % Initialize default SampEn method 

% Check data length, if < 34000 use Fast Implementation (introduced GDP) 
if length(data) < maxVecSize
     SampEnType = 'Fast'; 
end


% Loop through each window

% Loop through each timescale
% Note: i_tau == 1 is the original time series
for i_tau = 1:maxTau
    
    scaledData = coarsegrain(data,i_tau,cg_moment,cg_method); % Changed by GDP
    
    if ~constant_r
        scaledData = zscore(scaledData);
    end

    switch SampEnType
        case 'Fast'
           mse(i_tau) = fastSampen(scaledData, m, r);
        otherwise
           mse(i_tau) = sampenMaxim(scaledData, m, r); 
    end

end % end for loop






% REFERENCES

% Costa et al. "Multiscale entropy analysis of complex physiologic time
% series." Physical review letters 89.6 (2002): 068102

% Valencia et al. "Refined multiscale entropy: Application to 24-h holter 
% recordings of heart period variability in healthy and aortic stenosis
% subjects." IEEE Transactions on Biomedical Engineering 56, no. 9 (2009).

