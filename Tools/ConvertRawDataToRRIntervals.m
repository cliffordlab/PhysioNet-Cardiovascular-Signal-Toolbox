function [t,rr,jqrs_ann,SQIjw, StartIdxSQIwindows_jw] = ConvertRawDataToRRIntervals(ECG_RawData ,HRVparams, subjectID)
%   [t,rr,jqrs_ann,sqijs, StartIdxSQIwindows_jw] = ConvertRawDataToRRIntervals(ECG_RawData ,HRVparams, subjectID)  
%
%	OVERVIEW:
%       Load raw signal perfom QRS detection & Signal Quality Index SQI
%       and extract RR intervals from ECG signal 
%       (single lead ECG signal)
%
%   INPUT:
%       ECG_RawData : vector containing the 'raw' ECG signal (in mV) 
%       HRVparam    : struct of settings for hrv_toolbox analysis
%       subjectID   : name that identify the analyzed signal 
%
%   OUTPUT:
%       rr                     :  Vector containing RR interval
%       t                      :  Time indices of the rr interval data (seconds)
%       SQIjs                  :  Signal Quality Index values comparing jqrs and wqrsm 
%       StartIdxSQIwindows_jw  :  Indexes of SQI windows
%
%   Written by Giulia Da Poian (giulia.dap@gmail.com) 
% 

if nargin < 3
    error('Wrong number of arguments in ConvertRawDataToRRIntervals')
end

% run_sqrs and wqrsm detectors require ECG in digital values
% thus ECG data input that is in physical units will multiply the value
% GainQrsDetect

GainQrsDetect = 2000; % Default value for gain (adu/physical unit) 

if size(ECG_RawData,1)<size(ECG_RawData,2)
    ECG_RawData = ECG_RawData';
end

ECG_RawData = ECG_RawData(:,1); % If more the one leads use only the first

% QRS Dection 1 - jqrs
jqrs_ann = jqrs(ECG_RawData,HRVparams);

% QRS Detection 2 - sqrs (need single channel of ECG in digital values)
sqrs = run_sqrs(ECG_RawData*GainQrsDetect,HRVparams,0);

% QRS Detection 3 - wqrs
wqrs = wqrsm(ECG_RawData*GainQrsDetect);

% QRS SQI analysis
[sqijs, StartIdxSQIwindows_js] = bsqi(jqrs_ann(:),sqrs(:),HRVparams);
[SQIjw, StartIdxSQIwindows_jw] = bsqi(jqrs_ann(:),wqrs(:),HRVparams);


% Translate anotations to rr intervals
rr = diff(jqrs_ann./HRVparams.Fs);
t = jqrs_ann(2:end)./HRVparams.Fs;

%%  Export Annotations as ATR files

% Create a Folder for Annotations
WriteAnnotationFolder = [HRVparams.writedata filesep 'Annotation'];
if ~exist(WriteAnnotationFolder, 'dir')
   mkdir(WriteAnnotationFolder)
   fprintf('Creating a new folder: "Annotation", folder is located in %s \n',[pwd filesep WriteAnnotationFolder]);
end
addpath(WriteAnnotationFolder)

AnnFile = strcat(WriteAnnotationFolder, filesep, subjectID);
% Header File
write_hea(AnnFile, HRVparams.Fs, length(ECG_RawData), 'jqrs', 1, 0,'mV')
% ECG QRS and SQI
write_ann(AnnFile, HRVparams,'jqrs',jqrs_ann);
write_ann(AnnFile, HRVparams,'sqrs',sqrs);
write_ann(AnnFile, HRVparams,'wqrs',wqrs);
write_ann(AnnFile, HRVparams,'sqijs',sqijs);
write_ann(AnnFile, HRVparams,'sqijw',SQIjw); 





