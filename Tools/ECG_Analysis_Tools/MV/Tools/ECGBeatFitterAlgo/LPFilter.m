function y = LPFilter(x,fc)
%
% y = LPFilter(x,fc),
% Second order zero-phase Lowpass filter
%
% inputs:
% x: vector or matrix of input data (channels x samples)
% fc: -3dB cut-off frequency normalized by the sampling frequency
%
% output:
% y: vector or matrix of filtered data (channels x samples)
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.


k = .7;     % cut-off value
alpha = (1-k*cos(2*pi*fc)-sqrt(2*k*(1-cos(2*pi*fc))-k^2*sin(2*pi*fc)^2))/(1-k);
y = zeros(size(x));
for i = 1:size(x,1)
    y(i,:) = filtfilt(1-alpha,[1 -alpha],x(i,:));
end
