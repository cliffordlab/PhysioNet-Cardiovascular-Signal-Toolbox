% twa_demo.m, The following script demonstrates how to use the TWA code in
% the MVtoolbox to estimate TWAs using the modified moving average (MMA)
% method coupled with a non parametric surrogate statistical reshuffling
% method to reduce false positive TWA detections. The algorithm is
% referenced as follows,
%
%   Nemati S, Abdala O, Monasterio V, Yim-Yeh S, Malhotra A, Clifford GD. 
%   A nonparametric surrogate-based test of significance for T-wave alternans detection. 
%   IEEE Trans Biomed Eng. 2011;58(5):1356–1364. doi:10.1109/TBME.2010.2047859
%
% The main function is ComputeTWA. The output is stored in the TWAResult
% structure with the following fields with the expected output listed as follows,
%
%   Field                   Expected values
%   TWAResult.HR            [80.32,80.21,79.79,80.65,80.65,80.43,80,80,80.21,80.32,80,79.47]
%   TWAResult.VAlt          [38.84,40.10,37.81,37.38,37.70,38.14,39.32,39.74,38.31,40.39,39.07,42.20]
%   TWAResult.VAlt_Sig      [38.84,40.10,37.81,37.38,37.70,38.14,39.32,39.74,38.31,40.39,39.07,42.20]
%   TWAResult.Noise_Median  [4.28,3.92,4.40,3.75,3.68,4.14,3.00,3.21,4.15,3.95,3.78,7.20]
%   TWAResult.Noise_95      [11.19,10.75,9.61,12.67,11.69,9.72,9.15,9.84,10.29,11.03,9.89,19.40]
%   TWAResult.VAltPt        [132,116,115,125,116,122,130,117,91,111,101,159]
%   TWAResult.successful    1
%
% REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Ismail Sadiq
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  

% add dependencies within cardiovascular toolbox
addpath(genpath('../Tools/'));
addpath(genpath('../../../../../PhysioNet-Cardiovascular-Signal-Toolbox-master/'));
% load data and annotations
sig = load('./testdata/twa/sample_ecg_twam');
siginfo = readheader('./testdata/twa/sample_ecg_twam.hea');
Fs = siginfo.freq;
ecg = (sig.val - siginfo.adczero)./siginfo.gain; % Adjust signal according to gain and dc offset
clear sig siginfo
if (size(ecg,2) > size(ecg,1))
   ecg = ecg'; 
end

% load annotations
ann.QRSon = read_ann('./testdata/twa/sample_ecg_twa','qrson')'; ann.Q = read_ann('./testdata/twa/sample_ecg_twa','q')'; ann.R = read_ann('./testdata/twa/sample_ecg_twa','r')'; ann.S = read_ann('./testdata/twa/sample_ecg_twa','s')'; ann.QRSoff = read_ann('./testdata/twa/sample_ecg_twa','qrsoff')'; ann.Toff = read_ann('./testdata/twa/sample_ecg_twa','toff')';

% TWA evaluation
TWAResult = ComputeTWA(ecg,ann,Fs);