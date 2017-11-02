function rr_int = ApplyResamplingToRR(t_win,rr_0,HRVparams)

% OVERVIEW: Resampling of RR intervals before to apply methods that require 
%           resampled data for frequency domain analysis


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
