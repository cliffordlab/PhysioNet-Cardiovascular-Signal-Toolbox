function rr_int = ApplyResamplingToRR(t_win,rr_0,HRVparams)
%
%   rr_int = ApplyResamplingToRR(t_win,rr_0,HRVparams)
%	
%   OVERVIEW: Resampling of RR intervals before to apply methods that  
%             require resampled data for frequency domain analysis
%
%   INPUT:
%       t_win        - Vector containing start time of each window
%                     or ECG/PPG waveform  
%       rr_0         - vector of rr intervals to be resampled
%       HRVparams   - struct of settings for hrv_toolbox analysis that can
%                     be obtained using InitializeHRVparams.m function 
%                     HRVparams = InitializeHRVparams();
%   OUTPUT 
%       rr_int      - resampled rr intervals
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       This script written by Giulia Da Poian
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

sf = HRVparams.freq.resampling_freq; % (Hz) resampling frequency
ti = t_win(1):1/sf:t_win(end);       % time values for interp.
interp_method = HRVparams.freq.resample_interp_method;

% Chose the resampling method to use
switch interp_method 
    case 'cub'
        rr_int = interp1(t_win,rr_0,ti','spline')'; % cubic spline interpolation
    case 'lin'
        rr_int = interp1(t_win,rr_0,ti','linear')'; %linear interpolation
    otherwise
        warning('using cubic spline method')
        rr_int = interp1(t_win,rr_0,ti','spline')'; % cubic spline interpolation
end
