function Index = SVM_AFdetection_withoutTrainingModel(test_feature,test_class)
%******************************************************
% $ This function is used for calculating the SVM-based AF detection
% indices using the trained 9-fold models from MIT AF database. The
% determination of the AF is based on the voting resuts from the 9-fold
% models
%
% $ Reference: Q. Li, C. Y. Liu, J. Oster and G. D. Clifford. Chapter:
% Signal processing and feature selection preprocessing for classification
% in noisy healthcare data. In book: Machine Learning for Healthcare
% Technologies, Edition: 1st, Publisher: IET, Editors: David A. Clifton,
% 2016.
%
% Roberta Collaca, Julian Oster, Clifford
%
% $ Variable declaration: Input:
%   test_feature: the feature metrix from AF detection, you need to run
%       AF_feature function first to obtain the AF features, each row is a
%       feature vector from an RR interval time series
%   test_class: the labels of AF for the tested RR ineterval time series,
%       N*1 metrix, 1 for AF and 0 for non-AF
% Output:
%   Index: AF detection indices output, with the defined order of indices
%   of TP FN FP TN Se Sp Acc PPV NPV J in turn.
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

if nargin<2
    test_class = ones(length(test_feature),1);
end

% if strcmp(computer,'GLNXA64')
%     addpath('/home/adriana/libsvm-320/matlab');
% elseif strcmp(computer, 'MACI64')
%     addpath('/Users/adri/Documents/MATLAB/libsvm-3.20/libsvm-3.20/matlab/');
% end


model_all = load('models.mat');
meanx_all = load('means.mat');
stdx_all = load('stds.mat');

select=[2 3 4 5 6 7 12 14];
test_feature=test_feature(:,select);

predict_all=[];
for k = 1:9
    model=model_all.model_all{k};
    meanx=meanx_all.meanx_all{k};
    stdx=stdx_all.stdx_all{k};
    xtest=bsxfun(@minus, test_feature, meanx);
    xtest=bsxfun(@rdivide, xtest,stdx);
    [predict_all(:,k), accuracy_test, ytest_out] = svmpredict(test_class, xtest, model, '-b 1');
end

vote=sum(predict_all');
vote(find(vote<5))=0;
vote(find(vote>=5))=1;
vote=vote';

TP=length(find(test_class==1 & vote==1));
FP=length(find(test_class==0 & vote==1));
FN=length(find(test_class==1 & vote==0));
TN=length(find(test_class==0 & vote==0));
Se=TP/(TP+FN);
Sp=TN/(TN+FP);
Acc=(TP+TN)/length(test_class);
PPV=TP/(TP+FP);
NPV=TN/(TN+FN);
J=Se+Sp-1;
% Index=[TP FN FP TN Se Sp Acc PPV NPV J];

% Edited by Adriana N Vest
% Changed output Index to include only AF prediction, not accuracy
Index = [vote];
end



