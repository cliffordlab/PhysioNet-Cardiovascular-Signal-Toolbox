function [noise] = GenerateNoise(noisetype,SignalPower,SNR,N,varargin)
%
% noise =  GenerateNoise(noisetype,SignalPower,SNR,N,Param1,Param2,...),
% ECG noise generator
%
% Usage:
%       WN =  GenerateNoise(0,SignalPower,SNR,N,seed);
%       CN =  GenerateNoise(1,SignalPower,SNR,N,fs,beta,seed);
%       MX =  GenerateNoise(2,SignalPower,SNR,N,fs,[w_bw,w_em,w_ma],seed);
%
% inputs:
% noisetype
%       0:     white noise (WN)
%       1:     colored noise (CN)
%       2:     mixture of real baseline wander, electrode movements, muscle artifacts (MX)
% 
% SignalPower: The desired signal power. set to mean(x.^2) for the data vector x
% SNR: The desired SNR
% N: Number of samples
% fs: Sampling frequency required for noisetype = 1,2
% beta: Noise coloring factor required for noisetype = 1. beta = 0 (white noise),
%       beta = 1 (pink noise), beta = 2 (brown noise or random walk)
% seed(optional): Random seed for the noise vector. For noisetype =
%       2 seed is the initial random starting point in the real recorded
%       noises
% [w_bw,w_em,w_ma](optional): The weighting factors of BW, EM, and MA noise (only for noisetype = 2).
%
% output:
% noise: Column vector of noise
%
% This script is an editted version of the NoiseGenerator matlab function part of 
% the toolbox referenced below,
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com
%
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
%
% Last modified by Ismail Sadiq on 10/28/2019 to include randomized real
% noise subset selection.
%

switch noisetype
    case 0     % white noise
        if (nargin==5),
            randn('seed',varargin{1});
        end
        NoisePower = SignalPower / 10^(SNR/10);
        noise = sqrt(NoisePower)*randn(N,1);
        
    case 1     % colored noise
        fs = varargin{1};
        beta = varargin{2};
        %if (nargin==7),
        %    randn('seed',varargin{3});
        %end
        NoisePower = SignalPower / 10^(SNR/10);
        noise = ColoredNoise(sqrt(NoisePower),N,fs,beta);
        
    case 2     % mixture of real baseline wander, electrode movements, muscle artifacts
        if (nargin==7),
            randn('seed',varargin{3});
        end
        fs = varargin{1};
        w = varargin{2};
        w_bw = w(1);       % weight of baseline wander noise in the generated noise
        w_em = w(2);       % weight of electrode movement noise in the generated noise
        w_ma = w(3);       % weight of muscle artifact noise in the generated noise
        NoisePower = SignalPower / 10^(SNR/10);
        
        bw = csvread('bw.csv'); % read real baseline wander noise
        bw = bw(:,2);   % read channel 2 of recording
        em = csvread('em.csv'); % read real electrode movement noise
        em = em(:,2);   % read channel 2 of recording
        ma = csvread('ma.csv'); % read real muscle artifact noise
        ma = ma(:,2);   % read channel 2 of recording
        
        
        fs_noise = 360; % sampling frequency of each noise recording is 360 Hz
        bw = resample(bw',fs,fs_noise); % resample to same frequency as ECG
        bw = (bw-mean(bw))/std(bw); % scale noise recording to zero mean unit variance
        em = resample(em',fs,fs_noise);
        em = (em-mean(em))/std(em);
        ma = resample(ma',fs,fs_noise);
        ma = (ma-mean(ma))/std(ma);
        
        % weight each noise source as specified
        artifact = (w_bw*bw + w_em*em + w_ma*ma)/(w_bw + w_em + w_ma);
        
        if (nargin==7)
            % randomly select noise segment to add to ECG.
            n0 = randi(length(artifact)-N+1);
        else
            n0 = 1;
        end
        artifact = artifact(n0:N+n0-1)';
        % scale noise to appropriate level of power.
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);
end

%noise = noise(:);
end

