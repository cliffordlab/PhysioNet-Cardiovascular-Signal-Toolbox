function [NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt, SDNN, NNiqr, RMSSD, ...
          pnn50,btsdet,avgsqi,fbeatw, windows_all] = EvalTimeDomainHRVstats(...
          NN,tNN,sqi,HRVparams,windows_all,flagged_beats)
% [NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt,SDNN,NNiqr,RMSSD, ...
% pnn50,btsdet,avgsqi,fbeatw, windows_all] = EvalTimeDomainHRVstats(...
% windows_all,NN,tNN,sqi,settings,flagged_beats
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
%                                Additional columns can be included with
%                                additional sqi at the same timestamps
%               HRVparams      : struct of settings for hrv_toolbox analysis
%               windows_all    :
%               flagged_beats  : 
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
%               fbeatw      :
%               windows_all :
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Main script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%% 

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
if nargin <6 || isempty(flagged_beats)
    flagged_beats = zeros(length(NN),1);      
end


% Set Defaults

windowlength = HRVparams.windowlength;
dataoutput = HRVparams.timedomain.dataoutput;
alpha = HRVparams.timedomain.alpha;
threshold1 = HRVparams.timedomain.threshold1;
threshold2 = HRVparams.timedomain.threshold2;


% Preallocate arrays before entering the loop
NNmean     = zeros(1,length(windows_all));
NNmedian   = zeros(1,length(windows_all));
NNmode     = zeros(1,length(windows_all));
NNvariance = zeros(1,length(windows_all));
NNskew     = zeros(1,length(windows_all));
NNkurt     = zeros(1,length(windows_all));
NNiqr      = zeros(1,length(windows_all));
SDNN       = zeros(1,length(windows_all));
RMSSD      = zeros(1,length(windows_all));
pnn50      = zeros(1,length(windows_all));
btsdet     = zeros(1,length(windows_all));
avgsqi     = zeros(1,length(windows_all));
fbeatw     = zeros(1,length(windows_all));


%% Analyze by Window

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
        fbeats_win = flagged_beats(idx_NN_in_win);

        % Analysis of SQI for the window
        lowqual_idx = find(sqi_win(:,2) < threshold1);

        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < threshold2

            NNmean(i_win) = mean(nn_win);
            NNmedian(i_win) = median(nn_win);
            NNmode(i_win) = mode(nn_win);
            NNvariance(i_win) = var(nn_win);
            NNskew(i_win) = skewness(nn_win);
            NNkurt(i_win) = kurtosis(nn_win);
            NNiqr(i_win) = iqr(nn_win);
            SDNN(i_win) = std(nn_win); % SDNN should only be done on longer data segments

            % RMSSD
            RMSSD(i_win) = runrmssd(nn_win);

            % pNN50
            pnn50(i_win) = pnna(t_win, alpha);

            btsdet(i_win) = length(nn_win);
            avgsqi(i_win) = mean(sqi_win(:,2));
            fbeatw(i_win) = sum(fbeats_win);
                
          
        else
            NNmean(i_win) = NaN;
            NNmedian(i_win) = NaN;
            NNmode(i_win) = NaN;
            NNvariance(i_win) = NaN;
            NNskew(i_win) = NaN;
            NNkurt(i_win) = NaN;
            NNiqr(i_win) = NaN;
            SDNN(i_win) = NaN; 
            RMSSD(i_win) = NaN;
            pnn50(i_win) = NaN;
            btsdet(i_win) = length(nn_win);
            avgsqi(i_win) = mean(sqi_win(:,2));
            fbeatw(i_win) = sum(fbeats_win);
        end % end of conditional statements that run if SQI is above threshold2
        
    else
        NNmean(i_win) = NaN;
        NNmedian(i_win) = NaN;
        NNmode(i_win) = NaN;
        NNvariance(i_win) = NaN;
        NNskew(i_win) = NaN;
        NNkurt(i_win) = NaN;
        NNiqr(i_win) = NaN;
        SDNN(i_win) = NaN; 
        RMSSD(i_win) = NaN;
        pnn50(i_win) = NaN;
        btsdet(i_win) = NaN;
        avgsqi(i_win) = NaN;
        fbeatw(i_win) = NaN;
    end % end check for sufficient data
    
    % Print Results to File If Option Selected
    if dataoutput == 1
        if i_win == 1
            cd([HRVparams.writedata]) 
            fileID = fopen('generalstats.txt','w');
            formatSpec = '%6s %6s %6s %6s %6s %6s %6s %6s %7s %8s\n';
            fprintf(fileID,formatSpec,'NNmean','NNmedian','NNmode', ...
            'NNvariance','NNskew', 'NNkurt', 'NNiqr', 'SDNN', 'RMSSD', 'pnn50');
        end
        
        formatSpec = '%6.4f %6.4f %6.4f %6.4f %8.4f %8.4f %8.4f %8.4f %6.4f %6.4f \n';
        fprintf(fileID,formatSpec,NNmean(i_win),NNmedian(i_win),NNmode(i_win), ...
            NNvariance(i_win),NNskew(i_win) , NNkurt(i_win), NNiqr(i_win), ...
            SDNN(i_win), RMSSD(i_win), pnn50(i_win));
        
        if i_win == i_win(end)
            fclose(fileID);    
        end
    end

end % end of loop through window

end % end of function


