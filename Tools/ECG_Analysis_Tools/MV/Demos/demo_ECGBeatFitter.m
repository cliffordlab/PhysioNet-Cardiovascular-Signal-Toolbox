%
% Overview: Test ECGBeatFitter program for mean ECG beat extraction and
% parameter optimization. The script does the following:
% 1. Reads the ECG data in the sample_ecg_datam file in the 'Physionet-Cardiovascular-Signal-Toolbox/Tools/ECG_Analysis_Tools/MV/testdata/' subfolder.
% 2. Estimates the average ECG beat for the first fifteen seconds. 
% 3. The ECGBeatFitter algorithm is used to estimate the parameters of
% Gaussian functions that accurately estimate the shape of the average ECG
% beat. The number of Gaussians is specified by selecting points along the
% average beat in the ECGBeatFitter GUI.
% 4. The output of the script are the parameters for the Gaussians stored in the following variables:
%   
%   ai: contains the amplitudes of the Gaussians.
%   bi: contains the standard deviations of the Gaussians.
%   tetai: contains the phase of each Gaussian function.
%
% 5. These parameters may be estimated for each of the x, y and z
% components of a VCG and substituted into the generate_resp_modulated_ecg
% function in the Physionet-Cardiovascular-Signal-Toolbox to generate an artificial 
% VCG with a morphology similar to VCG used for estimating the parameters. 
% The Dower transform may be applied to generate 12 lead ECG from the VCG.
%
%
% ORIGINAL SOURCE AND AUTHORS: 
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com
% editted by Ismail Sadiq on 11/29/2020.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

clc
clear
close all;

% Load test ecg data
addpath(genpath('../../../../../PhysioNet-Cardiovascular-Signal-Toolbox-master/')); % add all dependencies from cardiovascular signal toolbox
sig = load('./testdata/mvm/sample_ecg_datam');
siginfo = readheader('./testdata/mvm/sample_ecg_datam.hea');
fs = siginfo.freq;
ecg = (sig.val - siginfo.adczero)./siginfo.gain; data=ecg(1:15000);% Adjust signal according to gain and dc offset
clear siginfo sig ecg;

t = (0:length(data)-1)/fs;

f = 1;                                          % approximate R-peak frequency

bsline = LPFilter(data,.7/fs);                  % baseline wander removal (may be replaced by other approaches)
%bsline = BaseLineKF(data,.5/fs);                % baseline wander removal (may be replaced by other approaches)

data1 = data-bsline;

%//////////////////////////////////////////////////////////////////////////
% Making the data noisy
SNR = 20;
SignalPower = mean(data1.^2);
NoisePower = SignalPower / 10^(SNR/10);
x = data1 + sqrt(NoisePower)*randn(size(data1));
%//////////////////////////////////////////////////////////////////////////

peaks = PeakDetection(x,f/fs);                  % peak detection

[phase, phasepos] = PhaseCalculation(peaks);     % phase calculation

teta = 0;                                       % desired phase shift
pphase = PhaseShifting(phase,teta);             % phase shifting

bins = round(fs/3);                                     % number of phase bins
[ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1); % mean ECG extraction 

OptimalParams = ECGBeatFitter(ECGmean,ECGsd,meanphase);               % ECG beat fitter GUI

% display the optimal parameters
L = length(OptimalParams)/3;
ai = OptimalParams(1:L)
bi = OptimalParams(L+1:2*L)
tetai = OptimalParams(2*L+1:3*L)

