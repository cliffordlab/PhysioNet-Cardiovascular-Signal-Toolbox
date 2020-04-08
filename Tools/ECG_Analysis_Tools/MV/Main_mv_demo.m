%Main_mv_demo.m
% OVERVIEW, Main script, runs the eval_mvm_twa script which calls the ComputeMVM and ComputeTWA in the Tools/MVM and Tools/TWA subfolders respectively. 
% 
% The expected output for this example script is given below.
% The results for MVM analysis are stored in the MVMResult structure and TWA analysis are stored in the TWAResult structure.
% The fields of each structure with expected values are listed below.
%
%   Structure.Field                 Expected value
%   MVMResult.energyinband_array    6.60e-07
%   MVMResult.sqi_array             0.9975
%   MVMResult.heart_rate_est_arr    80.11
%
%   TWAResult.HR            [80.32,80.21,79.79,80.65,80.65,80.43,80,80,80.21,80.32,80,79.47]
%   TWAResult.VAlt          [39.84,40.10,37.81,37.38,37.70,38.14,39.32,39.74,38.31,40.39,39.07,42.20]
%   TWAResult.VAlt_Sig      [39.84,40.10,37.81,37.38,37.70,38.14,39.32,39.74,38.31,40.39,39.07,42.20]
%   TWAResult.Noise_Median  [4.28,3.92,4.40,3.75,3.68,4.14,3.00,3.21,4.15,3.95,3.78,7.20]
%   TWAResult.Noise_95      [11.19,10.75,9.61,12.67,11.69,9.72,9.15,9.84,10.29,11.03,9.89,19.40]
%   TWAResult.VAltPt        [132,116,115,125,116,122,130,117,91,111,101,159]
%   TWAResult.successful    1
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

% Add path to the readheader function
addpath(genpath('../../../../PhysioNet-Cardiovascular-Signal-Toolbox-master'))
% Load data
sig = load('./Demos/testdata/twa/sample_ecg_twam');
siginfo = readheader('./Demos/testdata/twa/sample_ecg_twam.hea');
fs = siginfo.freq;
ecg = (sig.val - siginfo.adczero)./siginfo.gain; % Adjust signal according to gain and dc offset
clear sig siginfo
if (size(ecg,2) > size(ecg,1))
   ecg = ecg';
end

% Load annotations
ann.QRSon = read_ann('./Demos/testdata/twa/sample_ecg_twa','qrson')'; ann.Q = read_ann('./Demos/testdata/twa/sample_ecg_twa','q')'; ann.R = read_ann('./Demos/testdata/twa/sample_ecg_twa','r')'; ann.S = read_ann('./Demos/testdata/twa/sample_ecg_twa','s')'; ann.QRSoff = read_ann('./Demos/testdata/twa/sample_ecg_twa','qrsoff')'; ann.Toff = read_ann('./Demos/testdata/twa/sample_ecg_twa','toff')';
% Compute MVM and TWA in sample ECG.
%ann.Q = [];
[MVMResult,TWAResult] = eval_mvm_twa(ecg,ann,fs);
