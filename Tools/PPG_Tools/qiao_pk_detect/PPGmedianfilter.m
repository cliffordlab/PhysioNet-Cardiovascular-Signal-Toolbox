function xfilt = PPGmedianfilter(x, freq, freqd)
% FilterForTWA.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% Simple median filtration to remove baseline wander
% 12-01-2017 Modified by Giulia Da Poian: freqs as input parameter 
% Renamed to PPGmedianfilter (was medianfilter) to avoid conflict with
% other medianfilter functions
%	
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

% freqd = 32; % ECG

if nargin < 3
    freqd = 125; % PPG
end

factor = floor(freq / freqd);
y = x(1:factor:length(x));

halfwindclassic = floor(1.3 * freqd / 2);
 m = zeros(1,length(y));
for i = 1 : length(y)
    clear ss
    pstartc = max((i - halfwindclassic), 1);
    pendc = min((i + halfwindclassic), length(y));
    ss = sort(y(pstartc : pendc));
    m(i) = mean(ss(floor(0.375 * length(ss)) : floor(0.625 * length(ss))));    
end;

res = zeros(1, length(x));
for i = 1:length(x)
    res(i) = m(floor((i-1) / factor) + 1);
end;

xfilt = x - res';

return;
