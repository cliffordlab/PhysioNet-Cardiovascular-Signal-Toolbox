function [ noisy_s0 ] = AddNoise(s0,SNR,fs)
%OVERVIEW AddNoise: Adds baseline wander (BW), EMG and electrode movement artifacts (EA) to the signal
% at a specified SNR.
%   Inputs:
%       s0 - signal
%       SNR - signal to noise ratio
%       fs - sampling frequency
%   Outputs:
%       noisy_s0 - s0 with additive BW, EMG and EA
%   
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Ismail Sadiq
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.

Ch = size(s0, 2);   % No. ofN channels
N = size(s0, 1);    % length of data
noisy_s0 = zeros(N, Ch);

for k = 1:Ch
    signal_power = mean(s0(:,k).^2);    % get avg. power in lead
    %noise = GenerateNoise(0,signal_power,SNR,N,fs); % white noise
    %noise = GenerateNoise(1,signal_power,SNR,N,fs,1.5,1000); % colored noise
    noise = GenerateNoise(2,signal_power,SNR,N,fs,[1,1,1],1000); % real noise
    noisy_s0(:,k) = s0(:,k) + noise;
end

end

