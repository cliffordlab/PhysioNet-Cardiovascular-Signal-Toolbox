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
% Las updated: 2017.19.12 (by Giulia Da Poian) vectorized
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

dim = N-m;

% Fast implementation

t1 = repmat(x2(:,1),[1 dim]);
t2 = toeplitz([x2(1,1); flipud(x2(2:end,1))], x2(:,1));
t3 = repmat(x2(:,2),[1 dim]);
t4 = toeplitz([x2(1,2); flipud(x2(2:end,2))], x2(:,2));
d2max = max(abs(t1-t2),abs(t3-t4));
d2max = cell2mat(arrayfun(@(x) circshift(d2max(x,:),[1 dim-x+1]),(1:dim)','un',0));
Min_numerator = numel(find(d2max(dim,:)<=r));


% Old implementation (slower)
% % for i=1:N-m
% %     Min_numerator=0;
% %     for j=1:N-m
% %         d2max(i,j)=max(abs(x2(i,:)-x2(j,:)));
% %         if  d2max(i,j)<=r
% %             Min_numerator=Min_numerator+1
% %         end
% %     end
% % end


if N<20
    r=r+1;
else
    r=r+round(1000/fs);
end