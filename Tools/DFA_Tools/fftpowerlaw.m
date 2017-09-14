function y = fftpowerlaw(beta, n, seed)
% fftpowerlaw calculates power law data with exponent beta of length n
% fftpowerlaw(beta, n, seed) uses the specified seed.
% Edge effects (periodicities) are removed by only using the central half. 
% Copyright (c) 2005 by Patrick E. McSharry (patrick@mcsharry.net)

if nargin == 1
   n = 4096; 
end

if nargin == 2
   rand('seed', sum(100*clock) );
   randn('seed', sum(100*clock) );
end


% generate a white noise time series x 
N = 2*n;
x = randn(N,1);

% calculate the fft
Y = fft(x);
%k = [1:N]';
k = [1:N/2+1]';
Yp = ((k/(2*N)).^(-beta/2)).*Y(k);

% force symmetry of the spectrum 
Yp = [Yp; flipud(Yp(2:N/2))];
 
% Inverse fft and take the central half 
%y = real(ifft(Yp));
y0 = real(ifft(Yp));
ind = round(n/2); 
y = y0(ind+1:ind+n);


