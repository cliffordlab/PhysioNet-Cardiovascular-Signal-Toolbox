function features = AF_features(RR,fs)
%******************************************************
% $ This function is for calculating the common AF features
%
% $ Reference: Q. Li, C. Y. Liu, J. Oster and G. D. Clifford. Chapter:
% Signal processing and feature selection preprocessing for classification
% in noisy healthcare data. In book: Machine Learning for Healthcare
% Technologies, Edition: 1st, Publisher: IET, Editors: David A. Clifton,
% 2016.
%
% $ Variable declaration: 
%   Input:
%   RR:     RR interval time series (RR interval uses the unit of sample
%           points), RR time series are expected with the number of heart
%           beats range between 12 and 60. 
%   fs:     ECG sample rate 
%   Output:
%   features: AF features
%
% $ Author:  Chengyu Liu (chengyu.liu@emory.edu)
%           Department of Biomedical Informatics, Emory University, US
% $Last updated:  2016.10.28
% 
% %   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%
[nr,nc]=size(RR);
if nr<nc
    RR=RR';
end

M=floor(length(RR)/3);
RR=RR(1:3*M);
N=3*M;
if N<12 || N>60
    fprintf('Please input a RR interval time series with beat number between 12 and 60 \n');
end
m=1;
r=0.1;
%% time-domain features
mRR=mean(RR/fs);
minRR=min(RR/fs);
maxRR=max(RR/fs);
medHR=median(fs*60./RR);
SDNN=std(RR/fs);
dRR=diff(RR)/fs;
k=0;
sm=0;
for t=1:length(dRR)
    sm=sm+dRR(t).^2;
    if abs(dRR(t))>0.05
        k=k+1;
    end
end
PNN50=k/length(dRR)*100;
RMSSD=sqrt(sm/length(dRR));

%% frequency-domain features
p=6; % pole setting,
mhrf=1/mRR;
[P,f]=pburg(RR,p,[ ],[mhrf]); % using burg method
P(1)=P(2);
space=f(2)-f(1);
LF = sum(P(find( (0.04 <= f) & (f <= 0.15) )))*space;
HF = sum(P(find( (0.15  < f) & (f <= 0.40) )))*space;
Ratio_LH = LF/HF;
LFn = LF/(LF+HF);
HFn = HF/(LF+HF);

%% non-linear features
COSEn2 = COSEn(RR*1000/fs,m,30,fs);
COFMEn2 = COFMEn(RR*1000/fs,m,30,fs);
MAD= comput_MAD(RR*1000/fs);
AFEv = comput_AFEv2(RR/fs);

%% output features
features=[mRR minRR maxRR medHR SDNN PNN50 RMSSD LFn HFn Ratio_LH COSEn2 COFMEn2 MAD AFEv];
