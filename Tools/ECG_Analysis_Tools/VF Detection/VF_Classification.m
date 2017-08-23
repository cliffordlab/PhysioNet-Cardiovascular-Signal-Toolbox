function [VF_output,windows] = VF_Classification(tm,signal, s)
% VF_Classification
% Read WFDB format data file and classify VF / Non_VF
% input: 
%    record: record name of WFDB format data file
%
% output:
%    VF_output: VF (1) / Non_VF (0)
%
% This function need support of The WFDB Toolbox for MATLAB which can be
% download at http://physionet.org/physiotools/matlab/wfdb-app-matlab/
% to read WFDB format data
%
% This function need support of LIBSVM V3.11 for classification which can
% be download at http://www.csie.ntu.edu.tw/~cjlin/libsvm/
%
% Usage: VF_output = VF_Classification('cu01')
%
% Copyright (C) 2013 Qiao Li
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2013
% qiaolibme AT gmail DOT com
%
% Version 0.0.1
% Last updated: 06-27-2013
% 
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.
%
% Please cite this publication when referencing this program:
% Q Li, C Rajagopalan, GD Clifford, Ventricular fibrillation and tachycardia 
% classification using a machine learning approach, IEEE Transactions on 
% Biomedical Engineering, 61 (6), 1607-1613, 2014.
%
% April 17 2017
% Adriana N Vest Edits:
% Changed format of function to accept data directly from Matlab without
% using WFDB
%

% suppose the sampling frequency is 250
% if not, please resample the data to 250Hz using resample()
samp_freq = s.Fs;
powerline = 60;
window = 15; % 5 seconds analysis window
increment = 5;

% read data using WFDB_toolbox rdsamp
%[tm,signal,Fs] = rdsamp(record); % Commented out for Matlab only version

Fs = s.Fs;
tempdata = signal;
ECGleads = size(tempdata,2);
    
analysis_lead = 1; % analyze the first lead
analysis_data = tempdata(:,analysis_lead);
analysis_data(isnan(analysis_data)) = -5; % set NaN to -5 mV

% segments of data
n = floor(length(analysis_data)/samp_freq/window);
windows = 0:window:(n-1)*window;
load VF_classify_model

% 1 Hz high-pass to remove baseline wander
[b a] = butter(2,1.0/(samp_freq/2),'high');
% 30 Hz low-pass to remove high frequency noise
[b2 a2] = butter(2,30.0/(samp_freq/2));

for i=1:n
    begin_seg = (i-1) * samp_freq * window+1;
    end_seg = i * samp_freq * window;            
    data = analysis_data(begin_seg:end_seg);

    tempfilt = filtfilt(b,a,data);
    tempfilt = filtfilt(b2,a2,tempfilt);
    tempfilt = notchFilt(tempfilt,samp_freq,powerline);

    bandpass_data = band_pass(tempfilt);
    % calculate leakage
    leak = leakage(tempfilt);
    % calculate count2
    % (ANV added Samp Freq for calling counts.m)
    [count2] = counts(bandpass_data, s.Fs);
    % Classification by SVM, 1 for VF, 0 for Non_VF
    VF_output(i) = svmpredict(1, [leak,count2], model);
end
end

