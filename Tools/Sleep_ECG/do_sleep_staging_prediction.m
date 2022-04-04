function predict_label = do_sleep_staging_prediction(ECG,fs,class,trainingset) 
% input:
%   ECG: ECG waveform
%   fs: sampling frequency, 125Hz / 250Hz supported
%   class: categories to be classified
%       class = 1: four categories classification: Wake, REM sleep, NREM light sleep, NREM deep sleep
%       class = 2: three categories classification: Wake, REM sleep, NREM sleep (light + deep)
%       class = 3: three categories classification: Wake + REM sleep, NREM light sleep, NREM deep sleep
%       class = 4: two categories classification: Wake + REM sleep, NREM sleep (light + deep)
%   trainingset: the training data set used for the model, 'SHHSv1'(default) / 'SLPDB' / 'CinC2018tDB'
% 
% output:
% predict_label:
%   if OP.class = 1 then: predict_label = 1:Wake, 2:REM sleep, 3:NREM light sleep, 4:NREM deep sleep
%   if OP.class = 2 then: predict_label = 1:NREM sleep (light + deep), 2:REM sleep, 3: wake
%   if OP.class = 3 then: predict_label = 1:NREM deep sleep, 2:NREM light sleep, 3:Wake + REM sleep
%   if OP.class = 4 then: predict_label = 1:NREM sleep (light + deep), 2: Wake + REM sleep


if nargin<4
    trainingset='SHHSv1';
end

if nargin<3
    class=1;
end

if nargin<2
    fs=125;
end

OP.fs=fs;
if fs==125
    OP.lp_filter=@ecg_LP_filter_125_22;
    OP.hp_filter=@ecg_HP_filter_125_012;
elseif fs==250
    OP.lp_filter=@ecg_LP_filter_250_22;
    OP.hp_filter=@ecg_HP_filter_250_012;
else
    fprintf('Unsupported sampling frequency. Using "Filter Designer" toolbox to generate filters like "ecg_HP_filter_125_012.m" and "ecg_LP_filter_125_22.m"\n');
    return;
end
OP.class=class;
OP.trainingset=trainingset;
if strcmp(OP.trainingset,'SHHSv1')
    if class==1
        OP.epoch=39; % epoch of CNN used
        OP.model='CNN_SVM_model_SHHHSv1_class1.mat'; %'libsvm_noweight_shhs_HRV_kfold_output_resample_SVM_CNN_0_1_1.mat';
        OP.nclass=4; % number of categories
    elseif class==2
        OP.epoch=29;
        OP.model='CNN_SVM_model_SHHHSv1_class2.mat';
        OP.nclass=3;
    elseif class==3
        OP.epoch=47;
        OP.model='CNN_SVM_model_SHHHSv1_class3.mat';
        OP.nclass=3;
    elseif class==4
        OP.epoch=50;
        OP.model='CNN_SVM_model_SHHHSv1_class4.mat';
        OP.nclass=2;
    end
end

imdb = CRC_features_calculation(ECG,OP);

predict_label = CNN_SVM_prediction(imdb,OP);
