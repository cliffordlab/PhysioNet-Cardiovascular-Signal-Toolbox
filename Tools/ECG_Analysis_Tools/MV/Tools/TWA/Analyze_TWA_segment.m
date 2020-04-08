%function [TWARes.VAlt, TWARes.VAlt_Sig, TWARes.noise_median, TWARes.noise_95, TWARes.VAltPt, TWARes.threshold, HR] = Analyze_TWA_segment(ECG, Fs, ann)
function [TWARes] = Analyze_TWA_segment(ECG, Fs, ann)
% OVERVIEW, This function loads the parameters used for analysis in the Param structure, lowpass Filters the ecg and for power line interference if needed.
% It then calls the TWA_analyze script used for analyzing the ecg.
%
% INPUTS        MANDATORY       DESCRIPTION
%               ECG             N by M array of ECG data, with M channels
%                               and N datapoints in each channel.
%
%               Fs              Sampling frequency for the ecg, needs to
%                               1000 Hz.
%
%               ann             annotations for the ecg as structure
%                               variable. The variable contains the Q,R<S<T
%                               fiducial points.
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Ismail Sadiq
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  

% load twa analysis param
[~, Param] = set_twa_param();

% convert ecg to uV
ecg=ECG*1000;

% low pass filter ecg. also filter for power line interference if needed. 
ecg = Filter_TWA(ecg,60,40);

% perfrom TWA analysis
[TWARes] = TWA_analyze(ecg, ann, Param, 1000);
%TWARes.threshold = 0; %set theshold to 0;


end