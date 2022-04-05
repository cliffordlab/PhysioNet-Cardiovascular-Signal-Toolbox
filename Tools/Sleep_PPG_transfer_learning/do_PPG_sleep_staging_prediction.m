function [predict_label, predict_label_transfer, imdb] = do_PPG_sleep_staging_prediction(ppg,HRVparams,class,trainingset,SavePath,SigName) 
% input:
%           ppg : ppg waveform
%     HRVparams : HRVparams initialized
%         class : categories to be classified
%                 class = 1: four categories classification: (1) Wake, (2) REM sleep, (3) NREM light sleep, (4) NREM deep sleep
%                 class = 2: three categories classification: (1) NREM sleep (light+deep), (2) REM sleep, (3) Wake
%                 class = 3: three categories classification:  (1) NREM deep sleep, (2) NREM light sleep, (3) Wake + REM sleep,
%                 class = 4: two categories classification:  (1) NREM sleep (light + deep),  (2) Wake + REM sleep
%   trainingset : the training data set used for the model, 'SHHSv1'(default) / 'SLPDB' / 'CinC2018tDB'
%      SavePath :
%       SigName :
%
% output:
% predict_label:
%                if OP.class = 1 then: predict_label = 1:Wake, 2:REM sleep, 3:NREM light sleep, 4:NREM deep sleep
%                if OP.class = 2 then: predict_label = 1:NREM sleep (light + deep), 2:REM sleep, 3: wake
%                if OP.class = 3 then: predict_label = 1:NREM deep sleep, 2:NREM light sleep, 3:Wake + REM sleep
%                if OP.class = 4 then: predict_label = 1:NREM sleep (light + deep), 2: Wake + REM sleep
%
% 12-01-2018 Giulia : adding function for filtering ECG for a given Fs
% 12-02-2018 Giulia : adding imdb as additional output (to use for prediction in the future without 
%                     need to recompute the features, which is  time consuming)


if nargin<4
    trainingset='SHHSv1';
    SavePath='.';
    SigName='test';
end

if nargin<3
    class=1;
end

addpath('FeaturesExtraction')

OP.fs = HRVparams.Fs;

OP.class=class;
OP.trainingset=trainingset;
if strcmp(OP.trainingset,'SHHSv1')
    if class==1
        OP.epoch=39; % epoch of CNN used
        OP.model = ['CNN_SVM_model_SHHHSv1_class1.mat']; %'libsvm_noweight_shhs_HRV_kfold_output_resample_SVM_CNN_0_1_1.mat';
        OP.nclass=4; % number of categories
    elseif class==2
        OP.epoch = 29;
        OP.model = ['CNN_SVM_model_SHHHSv1_class2.mat'];
        OP.nclass = 3;
    elseif class==3
        OP.epoch = 47;
        OP.model = ['CNN_SVM_model_SHHHSv1_class3.mat'];
        OP.nclass=3;
    elseif class==4
        OP.epoch = 50;
        OP.model = ['CNN_SVM_model_SHHHSv1_class4.mat'];
        OP.nclass=2;
    end
end

if ~exist([SavePath filesep 'Features' filesep SigName '_imdb.mat'],'file')
    display('Features extraction. It may take a while. Please waiting...');
    imdb = Extract_PPG_feautures(ppg,HRVparams,SigName);
    if ~exist([SavePath filesep 'Features'], 'dir')
        mkdir([SavePath filesep 'Features']);
    end
    save([SavePath filesep 'Features' filesep SigName '_imdb.mat'],'imdb');
else
    load([SavePath filesep 'Features' filesep  SigName '_imdb.mat'],'imdb');
end

% prediction by ECG model
display('Prediction...');
predict_label = CNN_SVM_prediction(imdb,OP);

% transfer_learning_prediction

nclass=[4;3;3;2;2];
OP=[];
OP.fs = HRVparams.Fs;
OP.class=class;
OP.nclass=nclass(class);
balanced=1;

load('shhs_model_data_mean.mat');

load(['PPG_SVM_transfer_model_class' num2str(class) '_balanced_' num2str(balanced) '1']);
OP.model = ['PPG_SVM_transfer_model_class' num2str(class) '_balanced_' num2str(balanced) '1'];

[predict_label_transfer, predict_prob_transfer] = CNN_SVM_prediction_retrain(imdb,OP,stat.net,13,data_mean);


save([SavePath filesep SigName '_class' num2str(class) '.mat'],'predict_*');
