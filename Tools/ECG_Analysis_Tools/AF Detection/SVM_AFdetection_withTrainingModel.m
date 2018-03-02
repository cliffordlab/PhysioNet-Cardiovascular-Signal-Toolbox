function Index = SVM_AFdetection_withTrainingModel(training_feature,training_class,test_feature,test_class)
%******************************************************
% $ This function is used for calculating the SVM-based AF detection
% indices. You need to train the SVM model first using the training_feature
% and training_class, and then using the test_feature and test_class for calculating the AF
% detection incices
%
% $ Reference:
% Q. Li, C. Y. Liu, J. Oster and G. D. Clifford. Chapter: Signal processing
% and feature selection preprocessing for classification in noisy healthcare data. In book: Machine Learning for Healthcare Technologies, Edition: 1st, Publisher: IET, Editors: David A. Clifton, 2016.
%
% $ Variable declaration:
% Input:
% training_feature: the feature metrix from AF detection, you need to run
% AF_feature function first to obtain the AF features, each row is a
% feature vector from an RR interval time series
% training_class: the labels of AF for the trained RR ineterval time series, N*1
% metrix, 1 for AF and 0 for non-AF
% test_feature: the feature metrix from AF detection, you need to run
% AF_feature function first to obtain the AF features, each row is a
% feature vector from an RR interval time series
% test_class: the labels of AF for the tested RR ineterval time series, N*1
% metrix, 1 for AF and 0 for non-AF
% Output:
% Index: AF detection indices output, with the defined order of indices of TP FN FP TN Se Sp Acc PPV NPV J in turn. 
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information


addpath('E:\My algorithms\Matlab_SVM_tool\libsvm-3.20\matlab');

meanx=nanmean(training_feature);
stdx=nanstd(training_feature);

xtrain=bsxfun(@minus, training_feature, meanx);
xtrain=bsxfun(@rdivide, xtrain,stdx);

xtest=bsxfun(@minus, test_feature, meanx);
xtest=bsxfun(@rdivide, xtest,stdx);
model = svmtrain(training_class,xtrain, '-b 1');

[predict_all, accuracy_test, ytest_out] = svmpredict(test_class, xtest, model, '-b 1');

TP=length(find(test_class==1 & predict_all==1));
FP=length(find(test_class==0 & predict_all==1));
FN=length(find(test_class==1 & predict_all==0));
TN=length(find(test_class==0 & predict_all==0));
Se=TP/(TP+FN);
Sp=TN/(TN+FP);
Acc=(TP+TN)/length(test_class);
PPV=TP/(TP+FP);
NPV=TN/(TN+FN);
J=Se+Sp-1;
Index=[TP FN FP TN Se Sp Acc PPV NPV J];
end



