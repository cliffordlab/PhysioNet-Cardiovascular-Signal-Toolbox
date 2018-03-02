function [Min_numerator,r]=ReGet_min_numerator2(x1,x2,N,m,r,fs)
%******************************************************
% $ This function is usded to decrease r value if the many matchs are achieved. 
% $ Used for both COSEn and COFMEn. 
%
% $ Variable declaration: 
% x1 and x2 are the vector at dimentions of m and m+1
% N is time series length
% m is embedding dimension (usually m=1)
% r is changed threshold value
% fs sample rate
%
% $ Author:  Chengyu Liu (bestlcy@sdu.edu.cn) 
%           Institute of Biomedical Engineering,
%           Shandong University
% $Last updated:  2015.10.11
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

r=r-1;
for i=1:N-m
    Min_numerator=0;
    for j=1:N-m
        d2max(i,j)=max(abs(x2(i,:)-x2(j,:)));
        if  d2max(i,j)<=r
            Min_numerator=Min_numerator+1;
        end
    end
end
