function [s0] = generate_resp_modulated_ecg(fs, duration, resprate, heartrate,lambdascalar,snr)
%%
% Test program for generating synthetic abnormal multichannel ECGs.
% The current example is for the TWA abnormality.
%
% INPUTS        MANDATORY           DESCRIPTION
%               fs (Hz)             Sampling frequency of the synthetic ecg
%                                   to be generated. The script
%                                   has been tested for fs = 1000
%                                   Hz. For other fs values the synthetic 
%                                   ecg maybe resampled. 
%
%               duration (min.)     The duration of the synthetic ecg
%                                   segment in minutes.
%
%               resprate (rpm)      Resperation rate for synthetic ecg in 
%                                   respirations per minute (rpm).
%
%               heartrate (bpm)     Heart rate (HR) for synthetic ecg in
%                                   beats per minute (bpm).
%
%               lambdascalar        The lambda scalar used to generate the
%                                   T wave with elevated amplitude.
%
%               snr (dB)            Signal to noise ratio for synthetic ecg
%                                   generated.
% OUTPUTS
%               s0 (mV)             Sythetic ecg generated. N by Ch array,
%                                   where N is the data points in each channel 
%                                   and Ch is the number of channels.
%
% This script has been editted from the code provided in the toolbox referenced
% below.
%
% Open Source ECG Toolbox, version 2.0, April 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni [1], Gari Clifford [2]
% [1] Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, Grenoble, France, reza.sameni AT gmail DOT com
% [2] M.I.T., 77 Massachusetts Avenue, Cambridge MA 02139, gari AT mit DOT edu

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% Last editted on 10/28/2019 by Ismail Sadiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------- Dower Transform ----------

% declare the name of the Dower Transform Matrix (H) you want to use for
% mapping the VCG to the ECG leads produced. If this file doesn't
% exist, then one will be generated.
Hfile = 'T.dat';
% if we already have a Dower Transform Matrix, load it
if exist(Hfile);
    H = load(Hfile); % and load it if we have it.
    [a b] = size(H);
    if a<b; H=H'; end; % dimensionality consistency check.
    NumCh = size(H,1);
else
    H=[]; % else make H empty and we'll calcualte it form 1st principles later.
end
%% ----------------------- Initialization ---------------------------------

% ---------- dynamic range (number of bits used to save the data) ---------

% ---------- Noise parameters ----------

% ---------- Output Type ----------
Output_Type = 'ECG'; % VCG or ECG
%Num_Runs = 1;
% ---------- Extra Flags ----------
HRRampFlag = false; % makes the heart rate to ramp up in the middle
RSAFlag = false;
%RSFlag = true;
% RSFlag = false;
% ---------- General parameters ----------

len = 60*duration; % lenth of ECG in seconds
N   = fs*len;       % # of signal samples
%randoffset = 6;   % number of beats that the offset floor varies by
% ---------- Electrode pair locations ----------
ElecPos = [-10 5 15 ; -10 11 24 ; -10 0 23 ; -5 7 -7 ; -5 -1 -5 ; -10 10 18; -10 0 15; -10 10 15];
ElecNeg = [-10 5 24 ; -35 10 25 ; -10 10 24 ; 0 0 0 ; 0 0 0 ; -35 10 18; -10 10 15; -10 10 24];
% ---------- Heart location ----------
heartlocation = [-25 7 20];     % heart location with respect to the navel (coordinate reference)
% ---------- Heart Rate Parameters ----------

% ratio of  low frequency (LF) to high frequency (HF) in RR tachogram
lfhf  = 1/2;  stdhr = 5;% std of heart rate

k = 1;                      % dipole attenuation parameter
R0 = Rotate3D(0,0,0);       % dipole rotation matrices (tetax,tetay,tetaz)
Lambda = eye(3);
teta0 = -pi/3;              % initial phase of the ECG
%% ---------------------- End Of Initialization ---------------------------

Num_Runs = 1;

for run_indx = 1:Num_Runs
    for hr = heartrate
        for SNR = snr
            for loop = resprate % breathing rate scenarios
                
                if (loop == 0)
                    RSFlag = 0;     % respiration free for first value
                else
                    RSFlag = 1;     % respiration for others
                end
                
                tic; % set timer
                fprintf('Generating artificial ECG; BR = %i, hr= %i, #channels = %i, SNR = %i\n',loop,hr,NumCh,SNR);
                
                ramp =  40;       % ramp by which HR increases
                if (RSAFlag)
                    rr = rrprocess(ceil(hr*len/60), hr, lfhf, stdhr,1,0.1,loop/60);% generate an RR interval time series
                else
                    rr = rrprocess(ceil(hr*len/60), hr, lfhf, stdhr);% generate an RR interval time series
                end
                % and add a linear trend to RR intervals (in hr space, not rr space) +/- ramp/2 bpm
                hrv = 60./rr;
                if (HRRampFlag==true)
                    tmp=[1:1:length(rr)];tmp=ramp*tmp/max(tmp); tmp=tmp-mean(tmp); tmp=tmp';
                    % now make trnd suddenly ramp up in the middle, using tanh + some noise
                    tmp=((tanh(tmp)+1)/2)*ramp; tmp = tmp+(max(tmp)*rand(length(tmp),1)/30);
                    hrv=hrv+tmp; % add ramp to rr intervals in units of bpm
                end
                % convert to seconds
                F = 60./hrv;
                
                %===================== Normal Beat model =============================
                % Normal Beat model
                % alphai: structure contaning the amplitudes of Gaussian functions used for
                %           modeling the x, y, and z coordinates of the cardiac dipole
                % bi: structure contaning the widths of Gaussian functions used for
                %           modeling the x, y, and z coordinates of the cardiac dipole
                % tetai: structure contaning the phase of Gaussian functions used for
                %       modeling the x, y, and z coordinates of the cardiac dipole
                               
                tetai.x  = [-1.09  -0.83   -0.19  -0.07  0  0.06     0.22    1.2 1.42 1.68 2.9];
                alphai.x = [0.03   .08    -0.13    .85  1.11    0.75     0.06   0.1  lambdascalar*0.17 lambdascalar*0.39 lambdascalar*.03];
                bi.x     = [0.0906    0.1057    0.0453  0.0378    0.0332    0.0302    0.0378    0.6040 0.3020  0.1812 .5];
                
                tetai.y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8 1.58 2.9];
                alphai.y = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08 .014];
                bi.y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3 .5];
                                
                tetai.z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55 2.8];
                alphai.z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35 -.035];
                bi.z     = [.03  .12  .04         .4    .045       .05    .8 .4 .2 .4];
                                
                if isempty(H)
                    % ---------- Dower-like transform calculation ----------
                    
                    for i = 1:NumCh,
                        for j = 1:3,
                            H(i,j) = k* ((ElecPos(i,j)-heartlocation(j))/sqrt(sum((ElecPos(i,:)-heartlocation).^2))^3 - (ElecNeg(i,j)-heartlocation(j))/sqrt(sum((ElecNeg(i,:)-heartlocation).^2))^3);
                        end
                    end
                end
                % ---------- generate dipole ----------
                [DIP teta, hrt_rt] = DipoleGeneratorAbnormal(N,fs,F,alphai,bi,tetai,teta0);
                                
                % Rotate the dipole
                VCG = R0*Lambda*[DIP.x ; DIP.y ; DIP.z];
                
                if (RSFlag)
                    % simulating respiration induced amplitude modulation
                    n_resp_cycle = len/60 * loop;
                    ki = zeros(n_resp_cycle,1);ke = zeros(n_resp_cycle,1);
                    lami=zeros(n_resp_cycle,1); lame=zeros(n_resp_cycle,1);
                    ki(1)=0.35*fs; ke(1)=0.6*fs;
                    
                    fr = loop/60; % resp freq (breath per second)
                    lami = 20*fr/fs*ones(size(lami));
                    lame = 15*fr/fs*ones(size(lame));
                    for p=2:n_resp_cycle
                        ki(p) = ki(p-1)+fs/fr;
                        ke(p) = ke(p-1)+fs/fr;
                    end
                    Psi = 9;%*pi/180; % degree
                    QX=zeros(1,N);
                    for n=1:N
                        QX(n) = Psi*sum(1./((1+exp(-lami.*(n-ki))).*(1+exp(lame.*(n-ke)))));
                        Ri = Rotate3D(QX(n)*pi/180,QX(n)*pi/180,QX(n)*pi/180);
                        VCG(:,n)= Ri*VCG(:,n);
                    end
                end
                % ---------- convert VCG to ECG (clinical lead set) ----------
                
                if (strcmp(Output_Type,'ECG'))
                    s0 = (H*VCG)';
                else
                    s0 = (VCG)';
                    NumCh = 3;
                end
                % ---------- Scale the signal --------------
                for k=1:NumCh
                    max_amp(k) = max(abs(s0(:,k)));
                    s0(:,k) = (s0(:,k)/max_amp(k))*2;  % max at -+2 mV
                end
                
            end
        end
    end
end
end