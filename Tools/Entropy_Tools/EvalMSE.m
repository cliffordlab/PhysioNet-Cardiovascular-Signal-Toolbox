function mse = EvalMSE(NN,tNN,sqi,HRVparams,windows_all)

% mse = EvalMSE(NN,tNN,sqi,HRVparams,windows_all)
%
%   OVERVIEW:   This function returns MultiScale Entropy calculated on
%               input NN intervals for each window.
%
%   INPUT:      MANDATORY:
%               NN             : a single row of NN (normal normal) interval 
%                                data in seconds
%               tNN            : the time indices of the NN interval data 
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
%               mse            : vector of [max_tau, 1] doubles for each
%                                window
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


if isempty(HRVparams.MSE.windowlength)
    windowlength = length(NN);
else
    windowlength = HRVparams.MSE.windowlength*3600;
end

threshold1 = HRVparams.sqi.LowQualityThreshold;
threshold2 = HRVparams.RejectionThreshold;

m = HRVparams.MSE.patternLength;
r = HRVparams.MSE.RadiusOfSimilarity;
maxTau = HRVparams.MSE.maxCoarseGrainings;

cg_moment = HRVparams.MSE.moment;
cg_method = HRVparams.MSE.method;
constant_r = HRVparams.MSE.constant_r;

% Preallocate arrays (all NaN) before entering the loop
mse = nan(maxTau,length(windows_all));


%Analyze by Window

% Loop through each window of RR data
for i_win = 1:length(windows_all)
    % Check window for sufficient data
    if ~isnan(windows_all(i_win))
        % Isolate data in this window
        idx_NN_in_win = find(tNN >= windows_all(i_win) & tNN < windows_all(i_win) + windowlength);
        idx_sqi_win = find(sqi(:,1) >= windows_all(i_win) & sqi(:,1) < windows_all(i_win) + windowlength);
        
        sqi_win = sqi(idx_sqi_win,:);
        nn_win = NN(idx_NN_in_win);

        % Analysis of SQI for the window
        lowqual_idx = find(sqi_win(:,2) < threshold1);

        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < threshold2
            mse(:,i_win) = ComputeMultiscaleEntropy(nn_win, m, r, maxTau, [], cg_moment,cg_method,constant_r);
        end
        

    end % end check for sufficient data
    
end % end of loop through window

end % end of function


