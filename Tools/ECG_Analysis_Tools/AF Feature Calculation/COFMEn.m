function COFMEn = COFMEn(data, m, r,fs)
%******************************************************
% $ This function is usded for calcualting the cofficient of fuzzy measure entropy (COFMEn) for physiological signal time sequence. 
% $ Using COFMEn aims to improve the performance of COSEn. 
%
% $ Reference: Liu C Y, Li K, Zhao L N, Liu F, Zheng D C, Liu C C and Liu S
% T. Analysis of heart rate variability using fuzzy measure entropy.
% Computers in Biology and Medicine, 2013, 43(2): 100-108
%
% $ Variable declaration: 
% data is RR time series
% m is embedding dimension (usually m=1)
% r is threshold value (usually r=30 ms)
% local threshold r_l=0.2, global threshold r_g=0.2, 
% local weight of sequence segments' similarity n_l=2
% global weight of sequence segments' similarity n_g=2
%
% $ Author:  Chengyu Liu (bestlcy@sdu.edu.cn) 
%           Institute of Biomedical Engineering,
%           Shandong University
% $Last updated:  2015.10.15
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

% data=[108,108,109,109,108,108,110,107,109,109,109,109,108,109,108,108,108,108,107,108,108,108,108,107,109,107,108,107,108,108,108,108,108,108,108,109,108,109,108,108,109,108,108,108,109,107,108,107,107,108,108,108,109,107,108,108,108,108,108,108]*4;
% data=[109,108,108,109,109,109,108,109,108,108,109,108,108,109,108,108,108,109,108,109,109,109,108,108,109,108,108,109,109,109]*4;
% m=1;
% r=30;

N=length(data);
for i=1:N-m
    x1(i,:)=data(i:i+m-1);
    x2(i,:)=data(i:i+m);
end

ratio=0.4;
if N>20
%     Thr=ratio*N;
    Thr=5;
else
    Thr=5;
end


Min_numerator=0;
kk=0;
while Min_numerator<Thr
    [Min_numerator,r]=ReGet_min_numerator(x1,x2,N,m,r,fs);
    kk=kk+1;    
end
if kk==1 || kk>20
   if N<20
        r=r-1;
    else
        r=r-1000/fs;
    end
end

if Min_numerator==N-m
    while Min_numerator>=Thr
        [Min_numerator,r]=ReGet_min_numerator2(x1,x2,N,m,r,fs);
    end
    if r<0
        r=r+2;
    else 
        r=r+1;
    end
end

x=data;
r_l   = r;
r_g   = r;
n_l   = 2;
n_g   = 2;
%% m=2
D_l   = zeros(N-m,1); % initiate the mean local distance vector
D_g   = zeros(N-m,1); % initiate the mean global distance vector
for i = 1:N-m
    distance_l = zeros(N-m,1); % initiate the local distance vector for the ith vector
    distance_g = zeros(N-m,1); % initiate the global distance vector for the ith vector
    for j = 1:N-m
        if j==i
            d_l = 0;
            d_g = 0;
        else
        d_l = max(abs(x(i:i+m-1)-mean(x(i:i+m-1))-x(j:j+m-1)+mean(x(j:j+m-1))));
        d_g = max(abs(x(i:i+m-1)-x(j:j+m-1)));
        end
        distance_l(j) = exp(-(d_l.^n_l/r_l));
        distance_g(j) = exp(-(d_g.^n_g/r_g));
    end
    D_l(i) = sum(distance_l)/(N-m-1);
    D_g(i) = sum(distance_g)/(N-m-1);
end
Bm_l = mean(D_l);
Bm_g = mean(D_g);

%% m=m+1
m=m+1;
D_l   = zeros(N-m,1);
D_g   = zeros(N-m,1);
for i = 1:N-m
    distance_l = zeros(N-m,1);
    distance_g = zeros(N-m,1);
    for j = 1:N-m
        if j==i
            d_l = 0;
            d_g = 0;
        else
        d_l = max(abs(x(i:i+m-1)-mean(x(i:i+m-1))-x(j:j+m-1)+mean(x(j:j+m-1))));
        d_g = max(abs(x(i:i+m-1)-x(j:j+m-1)));
        end
        distance_l(j) = exp(-(d_l.^n_l/r_l));
        distance_g(j) = exp(-(d_g.^n_g/r_g));
    end
    D_l(i) = sum(distance_l)/(N-m-1);
    D_g(i) = sum(distance_g)/(N-m-1);
end
Am_l = mean(D_l);
Am_g = mean(D_g);

%% Calculate local and global fuzzy measure entropy
FuzzyLMEn = -log(Am_l/Bm_l);  % local fuzzy measure entropy
FuzzyGMEn = -log(Am_g/Bm_g);  % global fuzzy measure entropy
%% Generate fuzzy measure entropy
COFMEn  = FuzzyLMEn+FuzzyGMEn+2*log(2*r/1000)-2*log(mean(data)/1000);