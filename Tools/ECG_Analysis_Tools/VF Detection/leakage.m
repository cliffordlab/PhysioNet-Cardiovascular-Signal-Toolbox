function [leakage T] = leakage(data)
% Leakage calculation of ECG (VF filter)
% See: I. Jekova. Shock advisory tool: Detection of life-threatening 
% cardiac arrhythmias and shock success prediction by means of a common 
% parameter set.  Biomedical Signal Processing and Control, 2007(2):25-33
% See: Irena Jekova. Comparison of five algorithms for the detection of 
% ventricular fibrillation from the surface ECG. Physiol. Meas. 2000(21):
% 429-439
%
% Input: 
%  data: ECG data (one segment window)
%
% Output:
%  leakage: Leakage value  
%  T: period T
%
% Last modified by:
% Qiao Li, 4 Oct, 2012
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

data_n=length(data);
%PI=3.14159265359;

diff_d=diff(data);
t=2*pi*sum(abs(data))/sum(abs(diff_d));
t_2=round(t/2);

sum1=0;
sum2=0;
for i=t_2+1:data_n
    sum1=sum1+abs(data(i)+data(i-t_2));
    sum2=sum2+abs(data(i))+abs(data(i-t_2));
end

leakage=sum1/sum2;
T=t;
end