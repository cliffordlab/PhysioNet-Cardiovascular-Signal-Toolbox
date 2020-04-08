function [DIP, teta, hr]= DipoleGeneratorAbnormal(N,fs,rr,alphai,bi,tetai,teta0)
%
% [DIP teta]= DipoleGeneratorAbnormal(N,fs,rr,alphai,bi,tetai,teta0)
% Synthetic cardiac dipole generator using the 'differential form' of the
% dipole equations. Refer to references of the toolbox for further details.
%
% inputs:
% N: signal length
% fs: sampling rate
% rr: rr interval time series
% alphai: structure contaning the amplitudes of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% bi: structure contaning the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% tetai: structure contaning the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% teta0: initial phase of the synthetic dipole


%
% output:
% DIP: structure contaning the x, y, and z coordinates of the cardiac dipole
% teta: vector containing the dipole phase
%
%
% Open Source ECG Toolbox, version 2.0, April 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% Last update Shamim Nemati Feb 13 2009

rr = rr(:); t = [cumsum(rr)]; tp = linspace(0,sum(rr),N);
rrp = spline(t,rr,tp);
f=1./rrp; % make a heart rate (per second) time-series from rr intervals resampled to the current sampling frequency
hr = f*60;
w = 2*pi*f;     dt = 1/fs;
teta = zeros(1,N); X = zeros(1,N); Y = zeros(1,N); Z = zeros(1,N);

teta(1) = teta0;

for i = 1:N-1
    teta(i+1) = teta(i) + w(i)*dt; %dtheta/dt = w
    if(teta(i+1)>pi) % beat transition ------------------------------------
        teta(i+1) = teta(i+1) - 2*pi;
    end
    %------------------------------------------
    dtetaix = mod(teta(i) - tetai.x + pi , 2*pi) - pi;
    dtetaiy = mod(teta(i) - tetai.y + pi , 2*pi) - pi;
    dtetaiz = mod(teta(i) - tetai.z + pi , 2*pi) - pi;
    
    
    if(i==1),
        X(i) = sum(alphai.x .* exp(-dtetaix .^2 ./ (2*bi.x .^ 2)));
        Y(i) = sum(alphai.y .* exp(-dtetaiy .^2 ./ (2*bi.y .^ 2)));
        Z(i) = sum(alphai.z .* exp(-dtetaiz .^2 ./ (2*bi.z .^ 2)));
    end
    %------ Eq.(1) in Clifford et al. CinC2008 ------
    X(i+1) = X(i) - dt*sum(w(i)*alphai.x ./ (bi.x) .^ 2 .* dtetaix .* exp(-dtetaix .^2 ./ (2* bi.x .^ 2)));   % x state variable
    Y(i+1) = Y(i) - dt*sum(w(i)*alphai.y ./ (bi.y) .^ 2 .* dtetaiy .* exp(-dtetaiy .^2 ./ (2* bi.y .^ 2)));   % y state variable
    Z(i+1) = Z(i) - dt*sum(w(i)*alphai.z ./ (bi.z) .^ 2 .* dtetaiz .* exp(-dtetaiz .^2 ./ (2* bi.z .^ 2)));   % z state variable
    %------------------------------------------------
    
    
end
DIP.x = X;
DIP.y = Y;
DIP.z = Z;
