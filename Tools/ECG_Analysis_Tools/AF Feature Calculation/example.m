% **********************Zip file 2: AF feature calculation 
% AF feature
% calculation code package, which you can use to calculate the 14 AF
% features as the input feature metrix of SVM model. The AF 14 features
% inlcude the NFEn feature I sent you before.
% 
% ***********************

clear all
close all
clc

data=load('RR_example.txt');
features = AF_features(data(7,1:53),250);