function [alpha1, alpha2] = EvalDFA(NN,tNN,sqi,HRVparams,windows_all)

% [alpha1, alpha2] = EvalDFA(NN,tNN,sqi,HRVparams,windows_all)
%
%   OVERVIEW:   This function returns DFA scaling coefficients calculated on
%               input NN intervals for each window.
%
%   INPUT:      MANDATORY:
%               NN             : a single row of NN (normal normal) interval 
%                                data in seconds
%               tNN            : the time indices of the rr interval data 
%                                (seconds)
%               sqi            : (Optional )Signal Quality Index; Requires 
%                                a matrix with at least two columns. Column 
%                                1 should be timestamps of each sqi measure, 
%                                and Column 2 should be SQI on a scale from 0 to 1.
%               HRVparams      : struct of settings for hrv_toolbox analysis
%               windows_all    : vector containing the starting time of each
%                                windows (in seconds) 
%                
%   OUTPUT:     
%                  alpha1     : estimate of scaling exponent for short-term 
%                               fluctuations, minBoxSize <= n < midBoxSize
%                  alpha2     : estimate of scaling exponent for long-term 
%                              fluctuations, midBoxSize <= n <= maxBoxSize
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Giulia Da Poian    
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

% Verify input arguments
if nargin < 4
    error('Wrong number of input parameters');
end
if nargin < 5 || isempty(windows_all)
    windows_all = 0;   
end
if isempty(sqi) 
     sqi(:,1) = tNN;
     sqi(:,2) = ones(length(tNN),1);
end
% Set Defaults


if isempty(HRVparams.DFA.windowlength)
    windowlength = length(NN);
else
    windowlength = HRVparams.DFA.windowlength*3600;
end

SQI_LowQualityThresh = HRVparams.sqi.LowQualityThreshold;
RejectionThreshold = HRVparams.RejectionThreshold;

minBox = HRVparams.DFA.minBoxSize;
maxBox = HRVparams.DFA.maxBoxSize;
midBox = HRVparams.DFA.midBoxSize;

% Preallocate arrays (all NaN) before entering the loop
alpha1 = nan(length(windows_all),1);
alpha2 = nan(length(windows_all),1);

%Analyze by Window

% Loop through each window of RR data
for idxWin = 1:length(windows_all)
    % Check window for sufficient data
    if ~isnan(windows_all(idxWin))
        % Isolate data in this window
        idx_NN_in_win = find(tNN >= windows_all(idxWin) & tNN < windows_all(idxWin) + windowlength);
        idx_sqi_win = find(sqi(:,1) >= windows_all(idxWin) & sqi(:,1) < windows_all(idxWin) + windowlength);
        
        sqi_win = sqi(idx_sqi_win,:);
        nn_win = NN(idx_NN_in_win);

        % Analysis of SQI for the window
        lowqual_idx = find(sqi_win(:,2) < SQI_LowQualityThresh);

        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < RejectionThreshold
                alpha1(idxWin) = dfaScalingExponent(nn_win, minBox, midBox, 0);
                alpha2(idxWin) = dfaScalingExponent(nn_win, midBox, maxBox, 0);
        end
        

    end % end check for sufficient data
    
end % end of loop through window

end % end of function


