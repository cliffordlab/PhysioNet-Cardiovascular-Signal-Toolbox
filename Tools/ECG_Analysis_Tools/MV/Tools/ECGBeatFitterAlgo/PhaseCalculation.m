function [phase phasepos] = PhaseCalculation(peaks)
%
% [phase phasepos] = PhaseCalculation(peaks)
% ECG phase calculation from a given set of R-peaks.
%
% input:
% peaks: vector of R-peak pulse train
%
% outputs:
% phase: the calculated phases ranging from -pi to pi. The R-peaks are
% located at phase = 0.
% phasepos: the calculated phases ranging from 0 to 2*pi. The R-peaks are
% again located at phasepos = 0.
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.


phasepos = zeros(1,length(peaks));

I = find(peaks);
for i = 1:length(I)-1;
    m = I(i+1) - I(i);
    phasepos(I(i)+1:I(i+1)) = 2*pi/m : 2*pi/m : 2*pi;
end
m = I(2) - I(1);
L = length(phasepos(1:I(1)));
phasepos(1:I(1)) = 2*pi-(L-1)*2*pi/m:2*pi/m:2*pi;

m = I(end) - I(end-1);
L = length(phasepos(I(end)+1:end));
phasepos(I(end)+1:end) = 2*pi/m:2*pi/m:L*2*pi/m;

phasepos = mod(phasepos,2*pi);

phase = phasepos; 
I = find(phasepos>pi);
phase(I) = phasepos(I) - 2*pi;
