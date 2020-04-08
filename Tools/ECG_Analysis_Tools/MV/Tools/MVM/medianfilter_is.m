function [ signal_filter2 ] = medianfilter_is( signal, fs )
% [signal_filter2] = medianfilter_is( signal, fs )
%   OVERVIEW:   This function estimates the baseline wander signal in the ECG and returns the estimate in the signal_filter2 variable
%
%	INPUT: 	MANDATORY:
%               signal          : a single row of ECG data in samples.
%
%               fs              : sampling frequency for the ecg signal (Hz)
%
%
%   	OUTPUT:
%            	signal_filter2     : baseline wander estimate for the ecg in var signal.
%
%
%	REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Ismail Sadiq
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in the Documents 
%       folder of the Physionet-Cardiovascular-Signal-Toolbox.

orderfilter1 = floor(0.2 * fs);
orderfilter2 = floor(0.6 * fs);

signal_filter1 = medfilt1(signal, orderfilter1);
signal_filter2 = medfilt1(signal_filter1, orderfilter2);

end
