% mvm_demo

% OVERVIEW: This script provides a demo on how to execute the MVM function,
% part of the MVToolbox in the Physionet-Cardiovascular-Signal-Toolbox.
% This demo script does the following,
% 1. loads a sample ECG with annotations from the ./MVtoolbox/Demos/testdata subfolder
% 2. remove power line interference at 60 Hz and baseline wander signal
% from the ecg
% 3. Perform arrhythmia detection on each of the segment_size min analysis windows,
% avoid performing MVM analysis on windows with positive arrythmia
% detection
% 4. Compute MVM on segment_size min windows not affected by arrhythmia
% 5. store the output in the following variables,
%       variable            description (Units)
%       energyinband_array  QRS variability energy measured for each analysis window in the beatquency
%                           domain (beats^(-2))
%       sqi_array           The signal quality measured for each analysis
%                           window (unitless with value in the closed interval [0,1]. higher values mean cleaner signals)
%       heart_rate_estimate The estimate of the median heart rate for each
%                           analysis window is given in this array (beats per minute)
%
% The expected values for each of the output variables is as follows,
%       energyinband_array = [8.6362e-08,NaN,1.1978e-07]
%       sqi_array = [0.9988,NaN,1]
%       heart_rate_estimate = [80,NaN,79.893]
%
% REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Ismail Sadiq
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  

% load test ecg data
addpath('../../PhysioNet-Cardiovascular-Signal-Toolbox/Tools/ReadWrite_Files_Tools/'); % contains readheader function
sig = load('./testdata/sample_ecg_datam');
siginfo = readheader('./testdata/sample_ecg_datam.hea');
Fs = siginfo.freq;
ecg = (sig.val - siginfo.adczero)./siginfo.gain; % adjust signal according to gain and dc offset
clear siginfo sig;

% read QRSon, Q, R, S, QRSoff fiducial point annotations. Multiple
ann.QRSon = read_ann('./testdata/sample_ecgm','qrson'); ann.Q = read_ann('./testdata/sample_ecgm','q'); ann.R = read_ann('./testdata/sample_ecgm','r'); ann.S = read_ann('./testdata/sample_ecgm', 's'); ann.QRSoff = read_ann('./testdata/sample_ecgm', 'qrsoff');
%figure(1); plot(ecg); hold on; % check annotation if needed
%stem(ann.QRSon, ecg(ann.QRSon)); stem(ann.Q, ecg(ann.Q)); stem(ann.R, ecg(ann.R)); stem(ann.S, ecg(ann.S)); stem(ann.QRSoff, ecg(ann.QRSoff)); hold off

% read filter coefficients for power line interference (pli) filter
b = csvread('../Tools/pli_filter_coef_fs1000.csv');

% add path to helper function
addpath('../Tools/')
% add path to dependencies
addpath('../../PhysioNet-Cardiovascular-Signal-Toolbox/Tools/ECG_Analysis_Tools/PeakDetection_SQI/'); % path signal quality metric for ecg
addpath(genpath('../../PhysioNet-Cardiovascular-Signal-Toolbox/Tools/ECG_Analysis_Tools/PVCdetection/')); % path pvd detection tool

% remove baseline wander
baseline = medianfilter_is(ecg, Fs);
ecg_baselinefree = ecg - baseline;

% remove pli
% filter each signal for pli
ecg_filtered = filter(b,1,ecg_baselinefree); clear ecg_baselinefree baseline
% remove delay
delay = mean(grpdelay(b));
ecg_filtered(1:delay) = [];
%ecg_mv = ecg_filtered;

% detect arrhythmia
addpath(genpath('../../PhysioNet-Cardiovascular-Signal-Toolbox/Tools/ECG_Analysis_Tools/AF Feature Calculation/'));
addpath('../../PhysioNet-Cardiovascular-Signal-Toolbox/Tools/ECG_Analysis_Tools/AF Detection/');

% compute feature
minute_duration = 1*60*Fs;   % duration of 1 minute in samples
segment_size = 5;   % size of MVM analysis window in minutes
ii = 1; % index to keep track of each non overlapping segment_size minute window

energyinband_array = NaN(1,ceil(length(ecg_filtered)/(segment_size*minute_duration))); % array for storing the MVM energy evalauted for each segment_size minute analysis window, initialze to length of max possible non overlapping segment_size minute windows
sqi_array = energyinband_array; % array for storing the signal quality for each analysis window, same dimensions as energyinband_array
heart_rate_est = energyinband_array; % array for storing the heart rate estimate for each analysis window, same dimensions as energyinband_array
while (ii <= length(ecg_filtered))
    
    % use atrial fibirillation arrhythmia detection in 60 beat non overlapping windows for each segment_size min analysis window
    if (ii+segment_size*minute_duration-1 > length(ecg_filtered))
        % length of last analysis window is less than segment_size minutes, perform analysis till end of record 
        endi = length(ecg_filtered);
    else
        endi = ii+segment_size*minute_duration-1;
    end
    
    % get current 60 beat window for arrhythmia detection
    cur_win_r = ann.R(ann.R > ii); cur_win_r = cur_win_r(cur_win_r < endi);
    Index_test = NaN(1,ceil(length(cur_win_r)/60)); % variable for storing arrhythmia detection for each 60 beat window
    for jj = 1:60:length(cur_win_r)
        % get the current 60 beat window for arrhythmia detection
        if (jj+59 <= length(cur_win_r))
            cur_r = cur_win_r(jj:jj+59);
        else
            cur_r = cur_win_r(jj:end);
        end
        
        % if fewer than 12 beats or greater than 60 beats in the current
        % window do not perform arrhythmia detection
        if (length(diff(cur_r)) >= 12 && length(diff(cur_r)) <= 60)
            features = AF_features(diff(cur_r),Fs);
            Index_test(ceil(jj/60)) = SVM_AFdetection_withoutTrainingModel(features,1);
        else
            disp('Please input a RR interval time series with beat number between 12 and 60')
        end
        
    end
    
    if (nansum(Index_test) == 0)
        % if no arrhythmia detection in current segment_size min window, compute MVM
        
        % get qrson, q, r, s, qrsoff annotations for current window
        ann_win.R = cur_win_r-ii+1;
        ann_win.QRSon = ann.QRSon(ann.QRSon > ii); ann_win.QRSon = ann_win.QRSon(ann_win.QRSon < endi)-ii+1;
        ann_win.Q = ann.Q(ann.Q > ii); ann_win.Q = ann_win.Q(ann_win.Q < endi)-ii+1;
        ann_win.S = ann.S(ann.S > ii); ann_win.S = ann_win.S(ann_win.S < endi)-ii+1;
        ann_win.QRSoff = ann.QRSoff(ann.QRSoff > ii); ann_win.QRSoff = ann_win.QRSoff(ann_win.QRSoff < endi)-ii+1;
        
        % get ecg in current window
        ecg_filtered_win = ecg_filtered(ii:endi);
        % Morph Var QRS
        disp(['Evaluating QRS MVM for analysis window ' num2str(ceil(ii/(segment_size*minute_duration))) ' starting at time (s) ' num2str(ii/Fs)])
        debug = 0;
        normalize = 1;
        [energyinband_array(ceil(ii/(segment_size*minute_duration))),sqi_array(ceil(ii/(segment_size*minute_duration))),heart_rate_est(ceil(ii/(segment_size*minute_duration)))] = ComputeMVM(ecg_filtered_win,ann_win,Fs,segment_size,normalize);
        % analysis of current window successful, shift to next segment_size minute window
        ii = ii + segment_size*minute_duration;
    else
        % if arrhythmia detected in any consecutive 60 beats, shift analysis window by 1 minute. 
        disp('detected arrhythmia. shifting analysis window by 1 min.');
        ii = ii + minute_duration;
    end
        
    
end

