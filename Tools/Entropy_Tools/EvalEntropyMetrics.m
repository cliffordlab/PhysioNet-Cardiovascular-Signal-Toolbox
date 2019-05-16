function [SampEn, ApEn] = EvalEntropyMetrics(rr, t_rr, m, r, HRVparams, tWin, sqi)

% [SampEn, ApEn] = EvalEntropyMetrics(rr, t_rr, m, r, HRVparams, tWin, sqi)
%
% Overview
%	Calculates Sample Entropy and Approximate Entropy values of input data.
%
% Input
%   rr          - (seconds) rr intervals
%   t_rr        - (seconds) time of the rr intervals
%   m           - pattern length
%   r           - radius of similarity 
%   HRVparams   - struct of settings for hrv_toolbox analysis
%   sqi         - Signal Quality Index; Requires a matrix with at least two 
%                 columns. Column 1 should be timestamps of each sqi
%                 measure, and Column 2 should be SQI on a scale from 0 to 1
%   tWin        - vector containing the starting time of each
%                                windows (in seconds) 
% Output
%   SampEn      - 
%   ApEn        -
%
% Written by Giulia Da Poian <giulia.dap@gmail.com>
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%


if nargin<3
    Error('Not enough input arguments!');
end
if nargin<4 || isempty(HRVparams)
    windowlength = length(data);
    SQI_th = 0.9;          % SQI threshold
    WinQuality_th = 0.20;  % Low quality windows threshold
else
    windowlength = HRVparams.windowlength;
    SQI_th = HRVparams.sqi.LowQualityThreshold;        % SQI threshold
    WinQuality_th = HRVparams.RejectionThreshold;  % Low quality windows threshold
end
if nargin<5 || isempty(tWin)
    tWin = 0;
end
if nargin <6 || isempty(sqi)
    sqi(:,1) = t_rr;
    sqi(:,2) = ones(length(t_rr),1);
end

% Preallocation (all NaN)
SampEn = nan(length(tWin),1);
ApEn = nan(length(tWin),1);

% What Sample Entropy function?
SampEnType = 'Maxim'; % Initialize default SampEn method 
if windowlength < 34000; SampEnType = 'Fast'; end     
 
% Loop through each window of RR data
for iWin = 1:length(tWin)
    if ~isnan(tWin(iWin))
        % Isolate data in this window
        sqi_win = sqi( sqi(:,1) >= tWin(iWin) & sqi(:,1) < tWin(iWin) + windowlength,:);
        nn_win = rr( t_rr >= tWin(iWin) & t_rr < tWin(iWin) + windowlength );
        lowqual_idx = find(sqi_win(:,2) < SQI_th);  % Analysis of SQI for the window
        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < WinQuality_th
            
            nn_win = zscore(nn_win);  % normalization of the signal that replace the common 
                                      % practice of expressing the tolerance as r times the std
            % 1. Sample Entropy                                  
            switch SampEnType
                case 'Fast'
                   SampEn(iWin) = fastSampen(nn_win, m, r);
                otherwise
                   SampEn(iWin) = sampenMaxim(nn_win, m, r); 
            end
            % 2. Approximate Entropy
            ApEn(iWin) =  ApproxEntropy( nn_win, m, r);
        end 
    end 
end


