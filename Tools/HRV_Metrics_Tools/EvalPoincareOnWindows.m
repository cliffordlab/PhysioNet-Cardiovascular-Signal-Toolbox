function [SD1, SD2, SD1_SD2_ratio] = EvalPoincareOnWindows(rr, t_rr, HRVparams, tWin, sqi)

%   OVERVIEW:
%       Calculates SD1 SD2 and SD1_SD2_ratio features from Poincare plot
%       for each defined windows 
%   INPUT
%       rr           - (seconds) rr intervals
%       t_rr         - (seconds)  time stamp of rr intervals
%       HRVparams    - struct of settings for hrv_toolbox analysis
%       sqi          - Signal Quality Index; Requires a matrix with
%                      at least two columns. Column 1 should be
%                      timestamps of each sqi measure, and Column 2
%                      should be SQI on a scale from 0 to 1.
%       tWin         - Starting time of each windows to analyze 
%
%   OUTPUTS:
%       SD1           : (ms) standard  deviation  of  projection  of  the   
%                       PP  on  the line perpendicular to the line of 
%                       identity (y=-x)
%       SD2           : (ms) standard deviation of the projection of the PP  
%                       on the line of identity (y=x)
%       SD1_SD2_ratio : (ms) SD1/SD2 ratio
% 
% 
%   Written by: Giulia Da Poian <giulia.dap@gmail.com>
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox  
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

% Make vector a column
rr = rr(:);

if nargin < 4
    error('no data provided')
end
if nargin <5 || isempty(sqi)
    sqi(:,1) = t_rr;
    sqi(:,2) = ones(length(t_rr),1);
end

windowlength = HRVparams.windowlength;
SQI_th = HRVparams.sqi.LowQualityThreshold;        % SQI threshold
WinQuality_th = HRVparams.RejectionThreshold; % Low quality windows threshold

% Preallocation (all NaN)

SD1 = ones(length(tWin),1)*NaN;
SD2 = ones(length(tWin),1)*NaN;
SD1_SD2_ratio = ones(length(tWin),1)*NaN;

% Run PoincareMetrics by Windows
% Loop through each window of RR data
for i_win = 1:length(tWin)
    if ~isnan(tWin(i_win))
        % Isolate data in this window
        sqi_win = sqi( sqi(:,1) >= tWin(i_win) & sqi(:,1) < tWin(i_win) + windowlength,:);
        nn_win = rr( t_rr >= tWin(i_win) & t_rr < tWin(i_win) + windowlength );
        lowqual_idx = find(sqi_win(:,2) < SQI_th);         % Analysis of SQI for the window
        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < WinQuality_th
            [SD1(i_win), SD2(i_win), SD1_SD2_ratio(i_win)] = PoincareMetrics(nn_win);
        end % end of conditional statements run when SQI is adequate
    end % end of check for sufficient data
end % end of loop through windows



end % end of function
