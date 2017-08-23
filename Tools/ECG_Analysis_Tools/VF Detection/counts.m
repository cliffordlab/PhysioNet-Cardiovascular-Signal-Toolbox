function [count2] = counts(data,sampfreq) 
% See: Irena Jekova and Vessela Krasteva. Real time detection of 
% ventricular fibrillation and tachycardia. 2004 Physiol. Meas.
% 25:1167-1178
%
% Input: 
%  data: ECG data (one segment window)
%  sampfreq: sampling frequency
%
% Output:
%  count2: a metric for VF detection
%
% Last modified by:
% Qiao Li, 26 June, 2013
% qiaolibme AT gmail DOT com
% 
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.

if nargin<2
    sampfreq=250;
end

data_n = length(data);
data = abs(data);

count2=0;

for i=1:ceil(data_n/sampfreq)
    data_1s=data((i-1)*sampfreq+1:i*sampfreq);
    mean_d=mean(data_1s);
    count2=count2+length(find(data_1s>mean_d));
end
count2=count2/data_n;
end
