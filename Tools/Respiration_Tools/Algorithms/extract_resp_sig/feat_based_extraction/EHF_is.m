function [s_filt] = EHF_is()
%UNTITLED2 Summary of this function goes here


% ****Read file and resample



%   Detailed explanation goes here
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;
addpath('.\extract_resp_sig\filtering');

%% Load relevant data
if ~exist('data', 'var')
    %load([up.paths.data_load_folder, up.paths.data_load_filename]);
end

%% read mimic data
% mimicdb
[tm,signal,Fs,siginfo] = rdmat('221m');
% 221m - 2-6 min
% 224m - 
plot(tm, signal(:,5));
% resp signal 1-6 min
segdur = 1;
seglength = (segdur+5)*60;
tm_s = find(tm > 60); tm_s = tm_s(1);
tm_e = find(tm < seglength); tm_e = tm_e(end);
tm_seg = tm(tm_s:tm_e);
ecg = signal(tm_s:tm_e,1);
resp = signal(tm_s:tm_e,5);
figure(1); subplot(2,1,1); plot(tm_seg, ecg);
subplot(2,1,2); plot(tm_seg, resp);

% extract resp signal
ecg = resample(ecg,500,Fs); OP.fs = 500;
if (size(ecg,1) > size(ecg,2))
   ecg = ecg'; 
end

rel_data.v = ecg;
rel_data.fs = 500;
rel_data.method = 'real ECG';

%% Select appropriate filter characteristics
subj = 1;
curr_sig = 'ekg';
sig_type = 'ekg';
%eval(['rel_data = data(subj).' curr_sig ';']);
if strcmp(sig_type, 'ppg')
    filt_characteristics = up.paramSet.elim_vhf.ppg;
elseif strcmp(sig_type, 'ekg')
    %filt_characteristics = up.paramSet.elim_vhf.ekg;
    filt_characteristics.Fpass = 103.7;%up.paramSet.elim_vhf.ekg;
    filt_characteristics.Fstop = 98;
    filt_characteristics.Dpass = 0.05;
    filt_characteristics.Dstop = 0.01;
end

s = rel_data.v;
s_filt.fs = rel_data.fs;
s_filt.v = elim_vhfs(s, s_filt.fs, filt_characteristics);
s_filt.t = (1/s_filt.fs)*(1:length(s_filt.v));

%% eliminate mains frequencies
% only applicable to the ECG
filt_characteristics = [];
if strcmp(sig_type, 'ekg')
    filt_characteristics.fcuts = [41 47.2 52.8 59];%up.paramSet.elim_mains;
    filt_characteristics.Dpass = 0.05;
    filt_characteristics.Dstop = 0.01;
    s_filt.fs = rel_data.fs;
    s_filt.v = elim_mains(s_filt.v, s_filt.fs, filt_characteristics);
    s_filt.t = (1/s_filt.fs)*(1:length(s_filt.v));
end

% 'ekg_EHF = s_filt;'
end

