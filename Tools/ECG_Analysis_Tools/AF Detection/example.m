% **********************Zip file 1: SVM AF detection Two functions: 
% 1) SVM_AFdetection_withoutTrainingModel 
% 2) SVM_AFdetection_withTrainingModel
% 
% For function 1), you input the test AF feature metrix and label vector,
% and then obtain the evaluation indices for AF detection, such as true
% positive number, sensitivity, accuracy, etc. The AF detection uses the
% 9-fold SVM models trained with the MIT AF database. Three mat files
% included are the SVM parameters. 

% For function 2), you need input both the
% training AF feature metrix and label vector, as well as the test AF
% feature metrix and label vector, and then obtain the evaluation indices
% for AF detection. The function uses the training data to train the SVM
% and then get the model. For both 1) and 2), you need to call a famous
% open SVM classifier code libsvm-3.20, which is also included in the file.

clear all
close all
clc

training_feature = load('training_feature.txt');
training_class   = load('training_class.txt');
test_feature     = load('test_feature.txt');
test_class       = load('test_class.txt');

% Function 2
Index_train = SVM_AFdetection_withTrainingModel(training_feature, training_class, test_feature,test_class);

% Function 1
Index_test  = SVM_AFdetection_withoutTrainingModel(test_feature,test_class);


