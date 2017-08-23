function varargout = pwelchHRVtoolbox(x,varargin)
%PWELCH Power Spectral Density estimate via Welch's method.
%   Pxx = PWELCH(X) returns the Power Spectral Density (PSD) estimate, Pxx,
%   of a discrete-time signal, X, using Welch's averaged, modified
%   periodogram method.  When X is a vector, it is converted to a column
%   vector and treated as a single channel.  When X is a matrix, the PSD is
%   computed independently for each column and stored in the corresponding
%   column of Pxx.
%
%   By default, X is divided into the longest possible sections, to get as
%   close to but not exceeding 8 segments with 50% overlap. A modified
%   periodogram is computed for each segment using a Hamming window and all
%   the resulting periodograms are averaged to compute the final spectral
%   estimate. X is truncated if it cannot be divided into an integer number
%   of segments.
%
%   Pxx is the distribution of power per unit frequency. For real signals,
%   PWELCH returns the one-sided PSD by default; for complex signals, it
%   returns the two-sided PSD.  Note that a one-sided PSD contains the
%   total power of the input signal.
%
%   Note also that the default Hamming window has a 42.5 dB sidelobe
%   attenuation. This may mask spectral content below this value (relative
%   to the peak spectral content). Choosing different windows enables
%   you to make tradeoffs between resolution (e.g., using a rectangular
%   window) and sidelobe attenuation (e.g., using a Hann window). See
%   windowDesigner for more details.
%
%   Pxx = PWELCH(X,WINDOW), when WINDOW is a vector, divides each column of
%   X into overlapping sections of the same length as WINDOW, and then uses
%   the vector to window each section. If WINDOW is an integer, PWELCH
%   divides each column of X into sections of length WINDOW, and uses a
%   Hamming window of the same length. If the length of X is such that it
%   cannot be divided exactly into an integer number of sections with 50%
%   overlap, X is truncated. A Hamming window is used if WINDOW is omitted
%   or specified as empty.
%
%   Pxx = PWELCH(X,WINDOW,...,SPECTRUMTYPE) uses the window scaling
%   algorithm specified by SPECTRUMTYPE when computing the power spectrum.
%   SPECTRUMTYPE can be set to 'psd' or 'power':
%     'psd'   - returns the power spectral density
%     'power' - scales each estimate of the PSD by the equivalent noise
%               bandwidth of the window (in hertz).  Use this option to
%               obtain an estimate of the power at each frequency.
%   The default value for SPECTRUMTYPE is 'psd'.
%
%   Pxx = PWELCH(X,WINDOW,NOVERLAP) uses NOVERLAP samples of overlap from
%   section to section.  NOVERLAP must be an integer smaller than WINDOW if
%   WINDOW is an integer, or smaller than the length of WINDOW if WINDOW is
%   a vector. If NOVERLAP is omitted or specified as empty, it is set to
%   obtain a 50% overlap.
%
%   [Pxx,W] = PWELCH(X,WINDOW,NOVERLAP,NFFT) specifies the number of FFT
%   points used to calculate the PSD estimate.  For real X, Pxx has length
%   (NFFT/2+1) if NFFT is even, and (NFFT+1)/2 if NFFT is odd. For complex
%   X, Pxx always has length NFFT.  If NFFT is specified as empty, NFFT is
%   set to either 256 or the next power of two greater than the length of
%   each section of X, whichever is larger.
%
%   If NFFT is greater than the length of each section, the data is
%   zero-padded. If NFFT is less than the section length, the segment is
%   "wrapped" (using DATAWRAP) to make the length equal to NFFT. This
%   produces the correct FFT when NFFT is smaller than the section length.
%
%   W is the vector of normalized frequencies at which the PSD is
%   estimated.  W has units of radians/sample.  For real signals, W spans
%   the interval [0,pi] when NFFT is even and [0,pi) when NFFT is odd.  For
%   complex signals, W always spans the interval [0,2*pi).
%
%   [Pxx,W] = PWELCH(X,WINDOW,NOVERLAP,W) computes the two-sided PSD at the
%   normalized angular frequencies contained in the vector W.  W must have 
%   at least two elements.
%
%   [Pxx,F] = PWELCH(X,WINDOW,NOVERLAP,NFFT,Fs) returns a PSD computed as
%   a function of physical frequency.  Fs is the sampling frequency
%   specified in hertz.  If Fs is empty, it defaults to 1 Hz.
%
%   F is the vector of frequencies (in hertz) at which the PSD is
%   estimated.  For real signals, F spans the interval [0,Fs/2] when NFFT
%   is even and [0,Fs/2) when NFFT is odd.  For complex signals, F always
%   spans the interval [0,Fs).
%
%   [Pxx,F] = PWELCH(X,WINDOW,NOVERLAP,F,Fs) computes the two-sided PSD at 
%   the frequencies contained in the vector F.  F must have at least two
%   elements and be expressed in hertz.
%
%   [Pxx,F,Pxxc] = PWELCH(...,'ConfidenceLevel',P), where P is a scalar
%   between 0 and 1, returns the P*100% confidence interval for Pxx.  The
%   default value for P is .95. Confidence intervals are computed using a
%   chi-squared approach. Pxxc has twice as many columns as Pxx.
%   Odd-numbered columns contains the lower bounds of the confidence
%   intervals; even-numbered columns contain the upper bounds.  Thus,
%   Pxxc(M,2*N-1) is the lower bound and Pxxc(M,2*N) is the upper bound
%   corresponding to the estimate Pxx(M,N).
%
%   [...] = PWELCH(...,FREQRANGE)  returns the PSD over the specified range
%   of frequencies based upon the value of FREQRANGE:
%
%      'onesided' - returns the one-sided PSD of a real input signal X.
%         If NFFT is even, Pxx has length NFFT/2+1 and is computed over the
%         interval [0,pi].  If NFFT is odd, Pxx has length (NFFT+1)/2 and
%         is computed over the interval [0,pi). When Fs is specified, the
%         intervals become [0,Fs/2) and [0,Fs/2] for even and odd NFFT,
%         respectively.
%
%      'twosided' - returns the two-sided PSD for either real or complex
%         input X.  Pxx has length NFFT and is computed over the interval
%         [0,2*pi). When Fs is specified, the interval becomes [0,Fs).
%
%      'centered' - returns the centered two-sided PSD for either real or
%         complex X.  Pxx has length NFFT and is computed over the interval
%         (-pi, pi] for even length NFFT and (-pi, pi) for odd length NFFT.
%         When Fs is specified, the intervals become (-Fs/2, Fs/2] and
%         (-Fs/2, Fs/2) for even and odd NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after NOVERLAP.  The default value of FREQRANGE is 'onesided' when X
%      is real and 'twosided' when X is complex.
%
%   [...] = PWELCH(...,TRACE) uses TRACE to combine spectra computed
%   for each segment. TRACE can be set to 'mean', 'maxhold' or 'minhold':
%     'mean'    - Average spectrum over segments for each frequency bin
%     'maxhold' - Maximum spectrum over segments for each frequency bin
%     'minhold' - Minimum spectrum over segments for each frequency bin
%   The default value for TRACE is 'mean'.
%
%   PWELCH(...) with no output arguments plots the PSD estimate (in
%   decibels per unit frequency) in the current figure window.
%
%   EXAMPLE:
%      % Compute Welch's estimate of the periodogram of a 200 Hz sinusoid
%      % in noise using default values for window, overlap and NFFT.
%      Fs = 1000;   t = 0:1/Fs:.296;
%      x = cos(2*pi*t*200)+randn(size(t));
%      pwelch(x,[],[],[],Fs,'twosided');
% 
%   See also PERIODOGRAM, PCOV, PMCOV, PBURG, PYULEAR, PEIG, PMTM, PMUSIC.

%   Copyright 1988-2014 The MathWorks, Inc.

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice-Hall, 1997, pg. 15
%     [2] Monson Hayes, Statistical Digital Signal Processing and 
%         Modeling, John Wiley & Sons, 1996.

narginchk(1,9);
nargoutchk(0,3);

% look for psd, power, and ms window compensation flags
[esttype, varargin] = psdesttype({'psd','power','ms'},'psd',varargin);

% Possible outputs are:
%       Plot
%       Pxx
%       Pxx, freq
%       Pxx, freq, Pxxc
[varargout{1:nargout}] = welch(x,esttype,varargin{:});

% [EOF]
