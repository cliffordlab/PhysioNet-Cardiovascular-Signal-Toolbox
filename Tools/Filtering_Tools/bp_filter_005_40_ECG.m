function filtECG = bp_filter_005_40_ECG(ecg, fs)

% filtECG = filterECG(data, fs)
% apply a zerophase FIR filter 
% Equiripple FIR LP & HP filters (cascaded): 70dB 0.05-40Hz 1dB ripple. 
%
%
% Giulia Da Poian : giulia.dap@gmail.com


% Generate the filters
lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 35, ...
    'StopbandFrequency', 45, 'PassbandRipple', 1, ...
    'StopbandAttenuation', 70, 'SampleRate', fs, 'DesignMethod', 'equiripple');

hpFilt = designfilt('highpassfir', 'StopbandFrequency', .00128, ...
    'PassbandFrequency', 1.28, 'StopbandAttenuation', 70,...
    'PassbandRipple', 1, 'SampleRate', fs,'DesignMethod', 'equiripple');


den = 1;
numLP = lpFilt.Coefficients;
numHP = hpFilt.Coefficients;

% Data should be zeromeaned before calling this function
ecg = ecg-mean(ecg);

% low pass filter
ecg_lp = filtfilt(numLP,den,ecg);

% Don't high pass filter it if you don't want to remove the baseline 
% fluctuations due to resp, BP? and electrode noise?
filtECG = filtfilt(numHP,den,ecg_lp);
