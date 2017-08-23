% Cardiorespiratory Coupling Index (CRC Index)
%   [C,F] = crc_index(rr, ramp, nfft, fs, window, noverlap)
%   CRC is the product of the normalized output of 'csd_amp.m' and
%   'csd_phase.m'. The input rr and ramp should be an evenly sampled time
%   series.
%   CRC Index is the maximum CRC between 1.5 Hz and 5 Hz. If CRC index is
%   0, then F will be set to 0.
%   Usage:
%       Input:
%           rr = evenly sampled rr intervals.
%           ramp = evenly sampled ECG derived respiration
%           nfft = length of FFT.
%           Fs = sampling frequency.
%           window = hammming, hannning ...etc.
%           nooverlap = number of overlap.
%       Output:
%           C = cardiorespiratory coupling Index.
%           F = frequency at maximum CRC between 1.5 and 5 Hz.
%
%   Example:
%       Assuming rr and ramp are evenly sampled time series.
%       [C, F] = crc_index(rr,ramp,2^9, 30 ,hanning(200),100);
%
%   MATLAB Version 7.0.1.24704 (R14) Service Pack 1

function [c,f] = crc_index(rr, ramp, nfft, fs, window, noverlap)


[Pxy_amp, Fxy] = csd_amp(detrend(rr),detrend(ramp),nfft,fs,window,noverlap);
[Pxy_phase, Fxy] = csd_phase(detrend(rr),detrend(ramp),nfft,fs,window,noverlap);

Pxy_phase = abs(Pxy_phase);
Pxy_amp = abs(Pxy_amp);% Cross Spectral Density (Phase)

% get max amplitude and normalize
max_amp = max(Pxy_amp);
max_phase = max(Pxy_phase);
z1 = (Pxy_amp./max_amp).^2;
z2 = (Pxy_phase./max_phase).^2;
Cxy = z1.*z2;

%make sure Cxy and Fxy are column vectors
Cxy = Cxy(:);
Fxy = Fxy(:);

Cxy_diff = diff([0;Cxy]);% modified 3/21/2005

try
    % threshold
    Cxy_pos = Cxy_diff > 0;% modified 3/21/2005
    
    % locate index to look for maximum
    left  = find(diff([0; Cxy_pos])==1); % remember to zero pad at start
    right = find(diff([Cxy_pos; 0])==-1);  % remember to zero pad at end

    for(i=1:length(left))
        index = left(i):right(i);
        [C(i) maxloc(i)] = max(Cxy(index));
        F(i) = Fxy(index(maxloc(i)));
    end
    
    index2 = find(F > .15 & F < .4);
    [c, loc_c] = max(C(index2));
    
     % find frequency where Cxy is max
    f = F(index2(loc_c));
    
    
    if length(c) == 0
        c = 0;
        f = 0;
    end
catch
    c = 0;
    f = 0;
end

