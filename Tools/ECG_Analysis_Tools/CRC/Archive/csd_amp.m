% Cross Spectral Density (Amplitude)
%   [Pxy, Pxyc, f] = csd_phase(x,y, nfft,Fs,window, noverlap)
%   This mfile taken directly from Matlab mfile csd.m and was modified.
%   This a Welch's averaged cross periodogram method.
%   Usage:
%       Input:
%           x and y = signal vectors.
%           nfft = length of FFT.
%           Fs = sampling frequency.
%           window = hammming, hannning ...etc.
%           nooverlap = number of overlap.
%       Output:
%           Pxy = magnitude squared of cross spectral.
%           Pxyc = confidence interval (optional).
%           f = vector of frequencies.
%
%   Example:
%       Assuming RRI and RR_AMP are evenly sampled time series.
%       [Pxy_phase, Fxy] = csd_amp(detrend(RRI),detrend(RR_AMP),2^9, 2000 ,hanning(200),100);
%
%   MATLAB Version 7.0.1.24704 (R14) Service Pack 1



function [Pxy, Pxyc, f] = csd_phase(x,y, nfft,Fs,window, noverlap)

% compute CSD
window = window(:);
n = length(x);		% Number of data points
nwind = length(window); % length of window
if n < nwind    % zero-pad x , y if length is less than the window length
    x(nwind)=0;
    y(nwind)=0;  
    n=nwind;
end
x = x(:);		% Make sure x is a column vector
y = y(:);		% Make sure y is a column vector
k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
					% (k = fix(n/nwind) for noverlap=0)
index = 1:nwind;
KMU = k*norm(window)^2;	% Normalizing scale factor ==> asymptotically unbiased
% KMU = k*sum(window)^2;% alt. Nrmlzng scale factor ==> peaks are about right

Spec = zeros(nfft,1);
for i=1:k
    xw = window.*x(index);
    yw = window.*y(index);
    index = index + (nwind - noverlap);
    Xx = fft(xw,nfft);
    Yy = fft(yw,nfft);
    Xy2 = Yy.*conj(Xx);
    Spec = Spec + abs(Xy2);
end

% Select first half
if ~any(any(imag([x y])~=0)),   % if x and y are not complex
    if rem(nfft,2),    % nfft odd
        select = [1:(nfft+1)/2];
    else
        select = [1:nfft/2+1];   % include DC AND Nyquist
    end
    Spec = Spec(select);
else
    select = 1:nfft;
end
freq_vector = (select - 1)'*Fs/nfft;

p = [];
% find confidence interval if needed
if (nargout == 3)|((nargout == 0)&~isempty(p)),
    if isempty(p),
        p = .95;    % default
    end
    confid = Spec*chi2conf(p,k)/KMU;
end

Spec = Spec*(1/KMU);

% set up output parameters
if (nargout == 3)
   Pxy = Spec;
   Pxyc = confid;
   f = freq_vector;
elseif (nargout == 2)
   Pxy = Spec;
   Pxyc = freq_vector;
elseif (nargout == 1)
   Pxy = Spec;
elseif (nargout == 0)
   if ~isempty(p)
       P = [Spec confid];
   else
       P = Spec;
   end

    newplot;
    plot(freq_vector,10*log10(abs(P))), grid on
    xlabel('Frequency'), ylabel('Cross Spectrum Magnitude (dB)');
end
