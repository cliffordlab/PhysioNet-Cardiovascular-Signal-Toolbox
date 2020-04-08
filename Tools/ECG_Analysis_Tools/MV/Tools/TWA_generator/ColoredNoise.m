function n = ColoredNoise(sd,Len,fs,beta);
%
%  N = ColoredNoise(SD,Len,fs,beta),
%  Colored noise generation by frequency domain filtering of white noise.
%
% inputs:
% sd: the standard deviation of the generated noise
% Len: the noise vector length
% fs: the sampling rate
% beta: the noise color by assuming that the noise spectrum is propotional
% with 1/f^beta. beta = 0 corresponds to white noise.
% 
% output: colored noise vector
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.

len = Len + ceil(Len/10);
len = len + (4 - mod(len,4));

halflen = ceil(len/2);

s = sd*randn(len,1);
S = fft(s,len).^2;

f = [0:len-1]'*fs/len;
S(2:halflen+1) = S(2:halflen+1)./abs(f(2:halflen+1).^beta);
S(halflen+2:len) = S(halflen+2:len)./abs(f(halflen:-1:2).^beta);

n = real(ifft(sqrt(S),len,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A highpass filter for testing!
f1 = 0.1;
S(1) = 0;
k = 2:ceil(f1*len/fs);
S(k,:) = 0; S(len-k+2,:) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = n(halflen-ceil(Len/2):halflen+floor(Len/2)-1);
n = sd*(n-mean(n))/std(n);
n = n';