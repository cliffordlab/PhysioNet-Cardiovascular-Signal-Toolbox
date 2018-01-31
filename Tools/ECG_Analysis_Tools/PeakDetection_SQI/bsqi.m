function [F1, StartIdxSQIwindows] = bsqi(ann1,ann2,HRVparams)
%	OVERVIEW:
%       Measure SQI of ECG signals by comparing two peak detection
%       annotation files.
%   INPUT:
%       ann1         =  first annotation file
%       ann2         =  second annotation file
%       HRVparams    =  settings file 
%   OUTPUT:
%       output = sqi of each window
%   DEPENDENCIES & LIBRARIES:
%       run_sqi.m
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/Physionet-HRV-toolbox-for-MATLAB/tree/master/Tools
%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2017 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%
%   01-30-2018 - Modified by Giulia Da Poian (GDP): initialize F1 vector 
%               (vector of NaN), removed unused variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    HRVparams = initialize_HRVparams;
end

windowlength = HRVparams.sqi.windowlength;
threshold = HRVparams.sqi.TimeThreshold;
margin = HRVparams.sqi.margin;
fs = HRVparams.Fs;

ann1 = ann1(:)./fs;
ann2 = ann2(:)./fs;

endtime = max([ann1(end), ann2(end)]);
time = (1/HRVparams.Fs):(1/HRVparams.Fs):endtime;


%% 1. Create Windows
StartIdxSQIwindows = CreateWindowRRintervals(time, [], HRVparams,'sqi');

% Initialize Vectors
F1 = nan(1,length(StartIdxSQIwindows));

%% 2. Calculate SQI for Window
for seg = 1:length(StartIdxSQIwindows)
    % Check window for sufficient data
    if ~isnan(StartIdxSQIwindows(seg))
        % Isolate data in this window
        try
            idx_ann1_in_win = find(ann1 >= StartIdxSQIwindows(seg) & ann1 < StartIdxSQIwindows(seg) + windowlength);
            %idx_ann2_in_win = find(ann1 >= windows_all(i_win) & ann1 < windows_all(i_win) + windowlength);
        
            % Normalize timing of annotation data to the windows
            a1 = ann1(idx_ann1_in_win) - StartIdxSQIwindows(seg);
            a2 = ann2 - StartIdxSQIwindows(seg);
        
            F1(seg) = run_sqi(a1,a2,threshold,margin,windowlength,fs); % GDP : remove unused variables Se,PPV,Nb
        catch
        end
    end
end
