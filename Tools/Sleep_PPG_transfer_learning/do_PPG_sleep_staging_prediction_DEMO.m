record='test';
load('HRVparams_sleepstaging_PPG.mat');
load('shhs_model_data_mean.mat');

% read PPG data from WFDB/Matlab format
[tm,signal,Fs,siginfo]=rdmat(record);
class=1; % Four-class classification

[predict_label, predict_label_transfer, imdb] = ...
    do_PPG_sleep_staging_prediction(signal,HRVparams,class,'SHHSv1','.',record);
