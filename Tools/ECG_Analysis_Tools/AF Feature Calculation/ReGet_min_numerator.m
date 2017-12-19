function [Min_numerator,r]=ReGet_min_numerator(x1,x2,N,m,r,fs)
%******************************************************
% $ This function is usded to increase r value if the 'at least 5 match' is not achieved. 
% $ Used for both COSEn and COFMEn. 
%
% $ Variable declaration: 
% x1 and x2 are the vector at dimentions of m and m+1
% N is time series length
% m is embedding dimension (usually m=1)
% r is changed threshold value
% fs sample rate
%
% $ Author: Chengyu Liu (bestlcy@sdu.edu.cn) 
%           Institute of Biomedical Engineering,
%           Shandong University
% $Last updated:  2015.10.10
% Las updated: 2017.19.12 (by Giulia Da Poian) d2max preallocation

d2max = zeros(N-m, N-m);

for i=1:N-m
    Min_numerator=0;
    for j=1:N-m
        d2max(i,j)=max(abs(x2(i,:)-x2(j,:)));
        if  d2max(i,j)<=r
            Min_numerator=Min_numerator+1;
        end
    end
end
if N<20
    r=r+1;
else
    r=r+round(1000/fs);
end