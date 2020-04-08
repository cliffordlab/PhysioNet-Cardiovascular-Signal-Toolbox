% demo_wavedet_3D_detector.m
% OVERVIEW This script demonstrates the usage of the wavedet_3D_detector function.
% The wavedet_3D_detector function detects QRS onset, Q, R, S, J fiducial points for ECG required to perform MV analysis for ECG.
%
% The inputs and outputs to the main wavedet_3D_detector function are listed below.
% INPUTS        MANDATORY           DESCRIPTION
%               sig_all (uV)            Single channel ECG in mV.
%
%               ext_anot                The reference R peaks generated using jsqi. Provided to assist in beat detection. 
%
%               heasig                  Structure variable which contains
%                                       the following fields.
%               .nsamp                  Length of ecg data in sig_all in
%                                       samples.
%               .freq                   Sampling frequency for data in
%                                       sig_all.
%               .nsig                   Number of ecg channels in sig_all. Set to 1.                 
%
%               messages                            Strucuture variable with following field
%               .setup.wavedet.QRS_detection_only   Peforms fiducial point detection for the above listed 
%                                                   points of interest when set to 1.
%                                       
%
% OUTPUTS (of interest)
%               position                Structure variable which contains
%                                       the following fields of interest. 
%               .QRSon                  QRS onset location for each beat in
%                                       sig_all variable.
%               .Q                      Q point locations for each beat in
%                                       the sig_all variable.
%               .R                      R point locations for each beat in
%                                       the sig_all variable.
%               .S                      S point locations for each beat in
%                                       the sig_all variable.
%               .QRSoff                 QRS offset locations for each beat
%                                       in the sig_all variable.
%               .Toff                   T offset locations for each beat in
%                                       the sig_all variable.
%
% The detected fiducial points may be plotted along with the ecg to verify
% the demo executed correctly.
%
%   REF:
%       The wavedet_3D_detector algorithm has been adapted from the
%       opensource ecgkit toolbox referenced below,
%       Demski AJ, Llamedo Soria M. ecg-kit a Matlab Toolbox for Cardiovascular Signal Processing.
%       Journal of Open Research Software. 2016;4(1):e8. DOI: http://doi.org/10.5334/jors.86
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

% add toolbox dependencies
addpath(genpath('../../../../../PhysioNet-Cardiovascular-Signal-Toolbox-master/'));
addpath('../Tools/Annotation_generator/');
% load ecg data
sig = load('./testdata/mvm/sample_ecg_datam');
siginfo = readheader('./testdata/mvm/sample_ecg_datam.hea');
Fs = siginfo.freq;
ecg = (sig.val - siginfo.adczero)./siginfo.gain; % Adjust signal according to gain and dc offset
clear siginfo sig;

% generate annotations
HRVparams = InitializeHRVparams('MVanalysis');
[qrs_pos,sign,en_thres] = jqrs(ecg,HRVparams);
ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(ecg);
wavedet_config.setup.wavedet.QRS_detection_only = 0;
[ann,~,~] = wavedet_3D_detector(ecg, qrs_pos', ECG_header, wavedet_config );
% improve fiducial point detection
[detection_flg,ann] = improvfiducials(ann, Fs, ecg);    

if(~detection_flg)
    disp('Warning: Unable to improve fiducial point detection.')
end

% Plotting annotations for verification
figure(2); plot(ecg); hold on;
stem(ann.QRSon, ecg(ann.QRSon)); stem(ann.Q, ecg(ann.Q));
stem(ann.R, ecg(ann.R)); stem(ann.S, ecg(ann.S));
stem(ann.QRSoff, ecg(ann.QRSoff)); stem(ann.Toff(~isnan(ann.Toff)), ecg(ann.Toff(~isnan(ann.Toff)))); hold off
legend('ECG','QRSon','Q','R','S','QRSoff','Toff')
