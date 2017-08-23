function [filtdata]  = smp50Hzfilt(data, Fs, Fnotch)

% [filtdata]  = smp50Hzfilt(data, Fs, Fnotch);
% apply an FIR notch filter (coeffs designed by Neil Townsend)
% at Fnotch Hz for a sampling frequency of Fs Hz
% Defaults: Fnotch = 50.0;  Fs = 256;
%
% gari AT mit DOT edu ... G.D. Clifford 2004
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

% Data should be zeromeaned before calling this function

if nargin < 3
Fnotch  = 50.0;  % Hz    
end

if nargin < 2
Fs      = 256;  % Hz    
end

denominator = 1;
numerator = [1 1 1];
% numerator(1) = 1;
% numerator(3) = 1;
numerator(2) = -2.0*cos(2*pi*Fnotch/Fs);

% = 1.0/(2.0 + data->ds_smooth_fir_mult[ds][1]);

% make transfer function have unit gain
norm = sum(numerator);
numerator = numerator/norm;

% Forward-Reverse zero phase filter data
filtdata = filtfilt(numerator,denominator,data);
