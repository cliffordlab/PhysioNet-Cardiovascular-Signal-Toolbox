% demo_gen_twa_ecg.m
% OVERVIEW This script demonstrates the usage of the gen_twa_ecg function.
% The gen_twa_ecg function generates synthetic ecg for a given heart rate in beats per minute,
% respiration rate in respirations per minute, TWA amplitude in uV, sampling frequency in Hz and duration.
%
% The inputs and outputs to the main gen_twa_ecg function are listed below.
% INPUTS        MANDATORY           DESCRIPTION
%               twa_amp (uV)            The T wave alternan ampltiude in uV.
%
%               fs (Hz)                 Sampling frequency of the synthetic ecg
%                                       which will be generated. The script
%                                       has been tested for fs = 1000
%                                       Hz. For other fs values the synthetic ecg maybe resampled. 
%
%               duration (min.)         The duration of the ecg to be
%                                       sythesized in minutes.
%
%               hr (bpm)                The heart rate (HR) in beats per minute
%                                       (bpm).
%
%               rr (rpm)                The respiration rate (RR) in
%                                       respirations per minute (rpm).
%
%               SNR (dB)                Signal to noise ratio for sythetic
%                                       ecg signal.
%
% OUTPUTS
%               twa_ecg_n (mV)          Synthetic ecg with additive noise.
%                                       The TWA amplitude is specified by
%                                       the twa_amp variable.
%
%               twa_ecg (mV)            Noise free synthetic ecg. The TWA amplitude 
%                                       is specified by the twa_amp variable.
%
% The twa_ecg variable contains noise free ecg with 16 uV TWAs. The variable may be
% plotted to observe the TWAs.
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
%

% addpath to helper functions
addpath('../Tools/TWA_generator/');

fs = 1000; % sampling frequency Hz
hr = 80; % beats per minute
rr = 0; % respirations per minute
duration = 5; % minutes
twa_amp = 16; % microvolts
SNR = 30; % signal to noise ratio

[twa_ecg_n, twa_ecg] = gen_twa_ecg(twa_amp,fs,duration,hr,rr,SNR); % call the twa generation function





