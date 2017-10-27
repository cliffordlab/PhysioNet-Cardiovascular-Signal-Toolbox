function [NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt, SDNN, NNiqr, RMSSD, ...
          pnn50,btsdet,avgsqi,fdflag] = EvalTimeDomainHRVstats(NN,tNN,sqi,HRVparams,windows_all)

% [NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt,SDNN,NNiqr,RMSSD, pnn50
% btsdet,avgsqi] = EvalTimeDomainHRVstats(windows_all,NN,tNN,sqi,settings)
%
%   OVERVIEW:   This function returns time domain HRV metrics calculated on
%               input NN intervals.
%
%   INPUT:      MANDATORY:
%               NN             : a single row of NN (normal normal) interval data in seconds
%               OPTIONAL:
%               tNN            : the time indices of the rr interval data (seconds)
%               sqi            : Signal Quality Index; Requires a matrix with
%                                at least two columns. Column 1 should be
%                                timestamps of each sqi measure, and Column 2
%                                should be SQI on a scale from 0 to 1.
%               HRVparams      : struct of settings for hrv_toolbox analysis
%                
%   OUTPUT:     
%               NNmean      :  
%               NNmedian    :
%               NNmode      :
%               NNvariace   :
%               NNskew      : skewness
%               NNkurt      : kurtosis
%               SDNN        : standard deviation
%               NNiqr       :
%               RMSSD       :
%               pnn50       : the fraction of consecutive beats that differ by
%                             more than a specified time.
%               btsdet      :
%               avgsqi      :
%               fdflag      : 2 - Not enough high SQI data
%                             3 - Not enough data in the window to analyze
%                             5 - Success
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Adriana N. Vest     
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%
%   10-25-2017 Modified by Giulia Da Poian, convert results in ms and 
%   removed code for writing results on file
%

% Verify input arguments

if nargin< 1
    error('no input argument!!!!')
end
if nargin<2 || isempty(tNN)
        tNN = cumsum(NN);
end
if nargin<3 || isempty(sqi) 
     sqi(:,1) = tNN;
     sqi(:,2) = ones(length(tNN),1);
end
if nargin<4 || isempty(HRVparams) 
    HRVparams = initialize_HRVparams('demo');
end
if nargin<5 || isempty(windows_all)
    windows_all = CreateWindowRRintervals(tNN, NN, HRVparams);
end


% Set Defaults

windowlength = HRVparams.windowlength;
alpha = HRVparams.timedomain.alpha;
threshold1 = HRVparams.timedomain.threshold1;
threshold2 = HRVparams.timedomain.threshold2;


% Preallocate arrays (all NaN) before entering the loop
NNmean     = nan(1,length(windows_all));
NNmedian   = nan(1,length(windows_all));
NNmode     = nan(1,length(windows_all));
NNvariance = nan(1,length(windows_all));
NNskew     = nan(1,length(windows_all));
NNkurt     = nan(1,length(windows_all));
NNiqr      = nan(1,length(windows_all));
SDNN       = nan(1,length(windows_all));
RMSSD      = nan(1,length(windows_all));
pnn50      = nan(1,length(windows_all));
btsdet     = nan(1,length(windows_all));
avgsqi     = nan(1,length(windows_all));
fdflag     = nan(1,length(windows_all));

%Analyze by Window

% Loop through each window of RR data
for i_win = 1:length(windows_all)
    % Check window for sufficient data
    if ~isnan(windows_all(i_win))
        % Isolate data in this window
        idx_NN_in_win = find(tNN >= windows_all(i_win) & tNN < windows_all(i_win) + windowlength);
        idx_sqi_win = find(sqi(:,1) >= windows_all(i_win) & sqi(:,1) < windows_all(i_win) + windowlength);
        
        sqi_win = sqi(idx_sqi_win,:);
        t_win = tNN(idx_NN_in_win);
        nn_win = NN(idx_NN_in_win);

        % Analysis of SQI for the window
        lowqual_idx = find(sqi_win(:,2) < threshold1);

        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < threshold2

            NNmean(i_win) = mean(nn_win) * 1000; % compute and convert to ms
            NNmedian(i_win) = median(nn_win)* 1000; % compute and convert to ms
            NNmode(i_win) = mode(nn_win)* 1000; % compute and convert to ms
            NNvariance(i_win) = var(nn_win)* 1000; % compute and convert to ms
            NNskew(i_win) = skewness(nn_win)* 1000; % compute and convert to ms
            NNkurt(i_win) = kurtosis(nn_win)* 1000; % compute and convert to ms
            NNiqr(i_win) = iqr(nn_win)* 1000; % compute and convert to ms
            SDNN(i_win) = std(nn_win)* 1000; % compute and convert to ms % SDNN should only be done on longer data segments

            % RMSSD
            RMSSD(i_win) = runrmssd(nn_win)* 1000; % compute and convert to ms

            % pNN50
            pnn50(i_win) = pnna(nn_win, alpha); % 

            btsdet(i_win) = length(nn_win);
            avgsqi(i_win) = mean(sqi_win(:,2));
            
            fdflag(i_win) = 5; % 5 : 'sucess';
          
        else
            % 2: low SQI
            fdflag(i_win) = 2; 
        end % end of conditional statements that run if SQI is above threshold2
        
    else
        % 3: Not enough data in the window to analyze;
        fdflag(i_win) = 3; 
    end % end check for sufficient data
    
end % end of loop through window

end % end of function


