function [t,rr,jqrs_ann,SQIjw, StartSQIwindows_jw] = ConvertRawDataToRRIntervals(ECG_RawData ,HRVparams, subjectID)
%   [t,rr,jqrs_ann,sqijs, StartIdxSQIwindows_jw] = ConvertRawDataToRRIntervals(ECG_RawData ,HRVparams, subjectID)  
%
%   OVERVIEW:
%       Load raw signal, perform QRS detection & Signal Quality Index (SQI) analysis,
%       and extract RR intervals from ECG signal (single lead ECG signal).
%
%   INPUT:
%       ECG_RawData : vector containing the 'raw' ECG signal (in mV).
%       HRVparams   : struct of settings for hrv_toolbox analysis.
%       subjectID   : name that identifies the analyzed signal.
%
%   OUTPUT:
%       rr                  : (seconds) Vector containing RR intervals.
%       t                   : (seconds) Time of the RR interval data.
%       SQIjw               : Signal Quality Index values comparing jqrs and wqrs.
%       StartSQIwindows_jw  : time of SQI windows.
%
%   DEPENDENCIES & LIBRARIES:
%       PhysioNet Cardiovascular Signal Toolbox.
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox.
%
%   REFERENCE: 
%       Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%       Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%
%   REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox.
%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Giulia Da Poian (giulia.dap@gmail.com).
%       Dependent scripts written by various authors (see functions for details).
% 
%   EDITED BY
%       Weitao Tang
% 
%   COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if nargin < 3
    error('Wrong number of arguments in ConvertRawDataToRRIntervals')
end

% Default gain value for converting ECG data to digital values (adu/physical unit)
GainQrsDetect = 2000;

% Ensure ECG_RawData is a column vector
if size(ECG_RawData,1) < size(ECG_RawData,2)
    ECG_RawData = ECG_RawData';
end

% Use only the first lead if multiple leads are present
ECG_RawData = ECG_RawData(:,1); 

% QRS Detection 1 - jqrs
jqrs_ann = run_qrsdet_by_seg(ECG_RawData, HRVparams);

% QRS Detection 2 - sqrs (requires single channel of ECG in digital values)
sqrs_ann = run_sqrs(ECG_RawData * GainQrsDetect, HRVparams, 0);

% Check if sqrs_ann is empty, skip if detection failed
if isempty(sqrs_ann)
    warning('sqrs detection failed, skipping sqrs and using jqrs and wqrs only.');

    % QRS Detection 3 - wqrs
    wqrs_ann = wqrsm_fast(ECG_RawData * GainQrsDetect, HRVparams.Fs);

    % SQI analysis using jqrs and wqrs
    [SQIjw, StartSQIwindows_jw] = bsqi(jqrs_ann(:), wqrs_ann(:), HRVparams);
    SQIjs = []; % Set to empty if sqrs is not used
    StartSQIwindows_js = []; % Set to empty if sqrs is not used
else
    % If sqrs_ann is not empty, continue as usual
    wqrs_ann = wqrsm_fast(ECG_RawData * GainQrsDetect, HRVparams.Fs);

    % SQI analysis using jqrs, sqrs, and wqrs
    [SQIjs, StartSQIwindows_js] = bsqi(jqrs_ann(:), sqrs_ann(:), HRVparams);
    [SQIjw, StartSQIwindows_jw] = bsqi(jqrs_ann(:), wqrs_ann(:), HRVparams);
end

% Translate annotations to RR intervals
rr = diff(jqrs_ann ./ HRVparams.Fs);
t = jqrs_ann(2:end) ./ HRVparams.Fs;

%% Export Annotations as ATR files

% Create a folder for saving annotations
WriteAnnotationFolder = [HRVparams.writedata filesep 'Annotation'];
if ~exist(WriteAnnotationFolder, 'dir')
   mkdir(WriteAnnotationFolder)
   fprintf('Creating a new folder: "Annotation", located in %s \n', WriteAnnotationFolder);
end
addpath(WriteAnnotationFolder)

% Save annotations to file
AnnFile = strcat(WriteAnnotationFolder, filesep, subjectID);
write_hea(AnnFile, HRVparams.Fs, length(ECG_RawData), 'jqrs', 1, 0, 'mV');
write_ann(AnnFile, HRVparams, 'jqrs', jqrs_ann);
write_ann(AnnFile, HRVparams, 'sqrs', sqrs_ann);
write_ann(AnnFile, HRVparams, 'wqrs', wqrs_ann);

fakeAnnType = repmat('S', [length(SQIjs), 1]);
write_ann(AnnFile, HRVparams, 'sqijs', StartSQIwindows_js .* HRVparams.Fs, fakeAnnType, round(SQIjs * 100));

fakeAnnType = repmat('S', [length(SQIjw), 1]);
write_ann(AnnFile, HRVparams, 'sqijw', StartSQIwindows_jw .* HRVparams.Fs, fakeAnnType, round(SQIjw * 100));
