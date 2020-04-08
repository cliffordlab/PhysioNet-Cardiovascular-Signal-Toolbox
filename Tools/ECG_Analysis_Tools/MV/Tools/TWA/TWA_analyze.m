function [TWARes] = TWA_analyze(ecg, ann, Param, freq)
%OVERVIEW, THis function performs the TWA analysis and returns the results
%in the TWARes struct.
%
% INPUTS        MANDATORY
%
%               ecg (uV)        N by M array of ecg data. M channels with N
%                               datapoints for each channel.
%
%               ann             structure contains the Q,R,S,T-offset fiducial
%                               points.
%
%               freq (Hz)       samplg frequency of ecg, needs to be 1000
%                               Hz.
%
% OUTPUTS
%
%               TWARes.HR (bpm)             The average heart rate for each window analyzed
%
%               TWARes.VAlt (uV)            The T wave alternan amplitude estimate for each window analyzed.
%                               
%               TWARes.VAlt_Sig (uV)        The T wave alternan amplitude estimate for each window if it is statisitically significant compared to noise threshold.
%
%               TWARes.Noise_Median (uV)    The median of the gamma distribution used to estimate the noise distribution.
%
%               TWARes.Noise_95 (uV)        The 95th percentile of the gamma distribution used to estimate the noise distribution.
%
%               TWARes.successful           Flag which if set to 1 indicates if the TWA analysis is successful.
%
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
%

Param.stAdjIntv = floor(40 * freq / 1000); % interval of fiducial point adjustment (radius): depends on possible QS errors and QT variation.

[~,q,r,s,qs,stend,ecg]=get_ann(ecg, ann, Param, freq);%, []);

% removing last QRS to avoid problems with ecg length
s = s(1:length(s) - 1);
q = q(1:length(s));

% approximmate length of st segment
stlen = ApproximateSTLen(s, stend);

[TWARes] = TWA_MMA(Param, freq, ecg, q,r,s, stlen);
disp('TWA: done');


return;

