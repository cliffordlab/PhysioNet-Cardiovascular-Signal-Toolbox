function [PSD,F] = lombEC57(tNN,NN,bNorm)
%
% [PSD,F] = lombEC57(tNN,NN);
%
%   OVERVIEW:   Calculates the Lomb-Scargle normalized periodogram values
%               "PSD" as for input vectors "tNN" (time) and "NN" (observations).
%				Time stamps, tNN and amplitudes "NN" must be the same length.
%
%               Equivalent to the following native Matlab function 
%               (for bNorm=0): [PSD,F] = plomb(NN,tNN);
%               (for bNorm=1): [PSD,F] = plomb(NN,tNN,'normalized');
%
%   INPUT:      tNN  : the time indices of the rr interval data (seconds)
%               NN   : a single row of NN (normal normal) interval data in seconds
%               bNorm : normalize output; 1 for normalize, 0 for not
%
%   OUTPUT:     PSD  :
%               F    : frequencies at which PSD should be estimated
%                      default will be determined based on the maximal
%                      frequency detected in RR intervals (1/minimal interval)
%
%	Author:     Saeed Babaeizadeh, Advanced Algorithm Research Center, Philips Healthcare
%               based on plomb() in Matlab
%               August 2021
%
%%

% NN = NN-mean(NN); % subtract mean
L = length(NN);
nfft = 2^nextpow2(L); % Next power of 2 from length of signal
Fs = L/(tNN(end)-tNN(1));
% Calculate the number of unique points
NumUniquePts = ceil((nfft+1)/2);
% This is an evenly spaced frequency vector with NumUniquePts points.
f = (0:NumUniquePts-1)*Fs/nfft;

[PSD,F] = fasper(NN,tNN,4,f);

if bNorm
    varNN = var(NN);
    for i=1:length(PSD)
        if varNN~=0            % check for divide by zero
            PSD(i)=PSD(i)/(2*varNN);
        else
            PSD(i)=inf;
        end
    end
else
    PSD = PSD./Fs;
end

end

%% --------------------------------------------------------------------------
% Fast Lomb-Scargle Algorithm
function [wk2,f] = fasper(x,t,ofac,f)
%   ofac: oversampling factor,
%   The use of OFAC to interpolate or smooth a spectrum resembles the
%   zero-padding technique for FFT-based methods. P is again returned at
%   round(FMAX/FMIN) frequency points, but the minimum frequency considered
%   in this case is 1/(OFAC*N*Ts).

ofac=real(ofac(1,1));

n = numel(t);
tmin = t(1);
T = t(end)-tmin;
Ts = T/(n-1);

if numel(f) == 1
    fmax = f;
    nout = round(0.5*ofac*n);
    f = cast((1:nout)'/(n*Ts*ofac), 'like', x);
    hifac = cast(fmax(1,1)/max(f), 'like', x);
else
    hifac = cast(1, 'like', x);
end

MACC = 4;

nfreq = 2^nextpow2(ofac*hifac*n*MACC);
fndim = 2*nfreq;

% Initialize
wk1 = zeros(fndim,1);
wk2 = zeros(fndim,1);

% Extrapolate
fac = fndim/(n*Ts*ofac);
nfac = [1 1 2 6 24 120 720 5040 40320 362880];
nden = nfac(MACC);

for j = 1:n
    ck  = 1 + mod((t(j)-tmin)*fac,fndim);
    ckk = 1 + mod(2*(ck-1),fndim);
    wk1 = spread(x(j),wk1,fndim,ck ,MACC,nden);
    wk2 = spread(1   ,wk2,fndim,ckk,MACC,nden);
end

% Take the Fast Fourier Transforms
nout = round(0.5*ofac*hifac*n);
f = (1:nout)'/(n*Ts*ofac);

wk1 = fft(wk1);
rwk1 = real(wk1(2:nout+1));
iwk1 = imag(wk1(2:nout+1));

wk2 = fft(wk2);
rwk2 = real(wk2(2:nout+1));
iwk2 = imag(wk2(2:nout+1));

% Compute the Lomb-Scargle value for each frequency
hypo  = abs(wk2(2:nout+1));
if any(hypo==0)
    wk2 = lombscargle(x(:),f,t);
else
    hypo  = abs(wk2(2:nout+1));
    hc2wt = 0.5*rwk2./hypo;
    hs2wt = 0.5*iwk2./hypo;
    cwt   = sqrt(0.5+hc2wt);
    swt   = sign(hs2wt).*(sqrt(0.5-hc2wt));
    den   = 0.5*n + hc2wt.*rwk2 + hs2wt.*iwk2;
    cterm = (cwt.*rwk1 + swt.*iwk1).^2./den;
    sterm = (cwt.*iwk1 - swt.*rwk1).^2./(n-den);
    wk2  = cterm+sterm;
    wk2  = wk2(:);
end

end

%% --------------------------------------------------------------------------
% Extrapolate the frequency
function yy = spread(y,yy,n,x,m,nden)

ndenl = nden;

if x == round(x)
    yy(x) = yy(x) + y;
else
    ilo = min([max([floor(x-0.5*m+1),1]),n-m+1]);
    ihi = ilo+m-1;
    fac = (x-ilo)*prod(x-(ilo+1:ihi));
    yy(ihi) = yy(ihi) + y*fac/(ndenl*(x-ihi));
    for j = ihi-1:-1:ilo
        ndenl = double((ndenl/(j+1-ilo))*(j-ihi));
        yy(j) = yy(j) + y*fac/(ndenl*(x-j));
    end
    
end

end

%% --------------------------------------------------------------------------
% Conventional Lomb-Scargle Algorithm
function P = lombscargle(x,f,t)

nf = length(f);

if isrow(t)
    t1 = t(1,:).';
else
    t1 = t;
end

P = zeros(nf,1);

for i=1:nf
    wt = 2*f(i)*t1;
    swt = sinpi(wt);
    cwt = cospi(wt);
    
    Ss2wt = 2*cwt.'*swt;
    Sc2wt = (cwt-swt).'*(cwt+swt);
    
    wtau = 0.5 * atan2(Ss2wt,Sc2wt);
    
    swtau = sin(wtau);
    cwtau = cos(wtau);
    
    swttau = swt*cwtau - cwt*swtau;
    cwttau = cwt*cwtau + swt*swtau;
    
    swttauN = swttau.'*swttau;
    cwttauN = cwttau.'*cwttau;
    
    if abs(swttauN) < eps
        swttauN = swttauN + eps;
    end
    
    if abs(cwttauN) < eps
        cwttauN = cwttauN + eps;
    end
    
    P(i) = ((x.'*cwttau)^2)/(cwttauN) + ((x.'*swttau)^2)/(swttauN);
    
end
end

