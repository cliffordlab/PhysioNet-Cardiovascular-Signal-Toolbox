function sleep_stage = ECG_sleep_staging_prediction(ECGdata,fs,class)
% sleep stage prediction from ECG signal by model trained from SHHSv1 DB
%
% input:
%   ECGdata: ECG data
%   fs: sampling frequency of ECG
%   class: sleep staging type: class = 1, four types classification
%                              class = 4, two types classification
%
% output:
%   sleep_stage:  = 1: Wake, 
%                   2: REM sleep, 
%                   3: NREM light sleep, 
%                   4: NREM deep sleep, when class = 1
%   sleep_stage:  = 1: NREM sleep, 
%                   2: Wake+REM, when class = 4


classesn=[4 3 3 2 2];
numClasses=classesn(class);

% best parameters
filterSize=6;
balanced=1;
selectn = 3;

% 70 (10, sdnn), 73(14,avgsqi), 78(20, lfhf), 80(ac),81(dc),85(28, sampEn)
selects={[1:86],[1:69 71 72 74:77 79 82:84 86],[51:86],[51:69 71 72 74:77 79 82:84 86]};
used=selects{selectn};

load(['SHHS_ECG_model_' num2str(filterSize) '_' num2str(class) '_' num2str(balanced) '_' num2str(selectn) '_for_prediction'],'C','ka','parameters','hyperparameters','features_mean','features_std');  

% calculate features

features=ECG_sleep_staging_features_extraction(ECGdata,fs);

if ~isempty(features)
    doTraining=false;

    dlXTest=features;
    
    for l=1:length(features_mean)
            dlXTest(l,:)=bsxfun(@minus, dlXTest(l,:),features_mean(l));
            dlXTest(l,:)=bsxfun(@rdivide, dlXTest(l,:),features_std(l));
    end
    
    
    dlXTest=dlarray(dlXTest(used,:));
%     dlYTest = dlarray(s.Ydata_all);
    
    dimension = size(dlXTest,1);
    
    dlXTest=reshape(dlXTest,dimension,1,[]);
    dlYPred = model(dlXTest,parameters,hyperparameters,doTraining);
    dlYPred = softmax(dlYPred,'DataFormat','CBT');

%     YPred = onehotdecode(dlYPred,classes,1);
%     YTest = onehotdecode(dlYTest,classes,1);
%     
%     if isempty(predictions2)
%         predictions2=YPred(:);
%     else
%         predictions2 = [predictions2(:); YPred(:)];
%     end
%     yy=YPred == YTest;
%     predCorr2 = [predCorr2(:); yy(:)];   
    
    
    [~,YPred]=max(dlYPred);
    YPred=extractdata(YPred);
    sleep_stage=YPred(:);
    
%     YTest=extractdata(dlYTest);
%     YTest=YTest(:);
%     YTest=YTest+1;
% 
%     numFilters
%     C = confusionmat(YPred,YTest,'order',[numClasses:-1:1])
%     ka = kappa(C)
else
    sleep_stage = [];
end
