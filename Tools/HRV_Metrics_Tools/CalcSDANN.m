function [SDANN, SDNNI] = CalcSDANN(windows_all, tNN, NN, HRVparams)
% [SDANN, SDNNI] = CalcSDANN(windows_all, tNN, NN, HRVparams)
%
%   OVERVIEW:  Calculates SDANN and SDNNindex for all segments  
%              of a specified constant length with no overlap.
%
%   INPUT:     windows_all :
%              NN          : a single column of NN (normal normal) interval
%                            data in seconds
%              tNN         : the time indices of the rr interval data (seconds)
%              HRVparams   : settings struct, including segment length - the
%                            length of the interval requested 
%                            (it is defaulted at 5minutes == 300 seconds)
%              NOTE: seglength must be specified in the same units as <times>
%
%   OUTPUT:     SDANN : standard deviation of the averages of values
%               SDNNI : mean of the standard deviations of all values
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Main script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)  
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%

if nargin<3
    error('not enough input arguments!');   
end

segmentlength = HRVparams.windowlength;

for i = 1:length(windows_all)
    if ~isnan(windows_all(i))
        idx = find(tNN >= windows_all(i) & tNN < windows_all(i) + segmentlength);
    
        nn_win = NN(idx);
    
        sm(i) = mean(nn_win);     % mean of each segment
        sstd(i) = std(nn_win);    % stdev of each segment
    else
        sm(i) = NaN;
        sstd(i) = NaN;
    end
end

% Remove NaN values from calculation
seg_mean = sm;
seg_stdev = sstd;
idxrem = find(isnan(sm));
seg_mean(idxrem) = [];
seg_stdev(idxrem) = [];
    
SDANN = std(seg_mean);              % stdev of means of each segment
SDNNI = mean(seg_stdev);            % mean of stdevs from each segment

end