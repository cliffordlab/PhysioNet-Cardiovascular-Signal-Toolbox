function avgsqi = WindowAvgSqi(sqi,tNN,HRVparams)

%   avgsqi = WindowAvgSqi(sqi,tNN,HRVparams)
%
%   OVERVIEW:   
%
%   INPUT:      MANDATORY:
%               sqi         : Signal Quality Index
%               tNN         : the time indices of the rr interval data (seconds)
%               HRVparams   : struct of settings for hrv_toolbox analysis
%                               
%   OUTPUT:     
%               avgsqi      : average value of SQIs
%               
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
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

[~,col] = size(sqi);

%% Create windows for analysis of all windows
windows_all = [];


for j = 1:(floor(tNN(end))-(HRVparams.windowlength))/HRVparams.overlap
    indicies = [];
    try
        indicies = find(tNN >= (j-1)*HRVparams.overlap & tNN < HRVparams.windowlength+(j-1)*HRVparams.overlap);
        % windows_all: first coloumn is start time of sliding window
        % and second column is end time of sliding window
        windows_all(j,1) = tNN(indicies(1));
        windows_all(j,2) = tNN(indicies(end));
    catch
        windows_all(j,1) = NaN;
        windows_all(j,2) = NaN;
    end
end

% Replace windows that aren't equal to variable windowlength (with a tolerance) 
% with NaN values 
meas_win_length = windows_all(:,2)-windows_all(:,1);
too_short = find(meas_win_length < (HRVparams.windowlength * (1-HRVparams.win_tol)));
for k = 1:length(too_short)
    windows_all(too_short(k),1) = NaN;
    windows_all(too_short(k),2) = NaN;
end

% Loop through each window of RR data
for i_win = 1:length(windows_all)
    % Check window for sufficient data
    if ~isnan(windows_all(i_win,1))
        %idx_NN_in_win = find(tNN >= windows_all(i_win,1) & tNN < windows_all(i_win,1) + s.windowlength);
        idx_sqi_win = find(sqi(:,1) >= windows_all(i_win,1) & sqi(:,1) < windows_all(i_win,1) + HRVparams.windowlength);

        % Isolate data in this window
        for k = 1:(col-1)
            sqi_win = sqi(idx_sqi_win,k+1);
            avgsqi(i_win,k) = mean(sqi_win);
        end
    end   
end
end % end of function