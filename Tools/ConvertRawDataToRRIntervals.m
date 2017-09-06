function [t,rr,jqrs_ann,sqijw] = ConvertRawDataToRRIntervals(ECG_RawData ,HRVparams, subjectID)
%   [t,rr,jqrs_ann,sqijw] = ConvertRawDataToRRIntervals(HRVparams)  
%
%	OVERVIEW:
%       Load raw signal perfom QRS detection & Signal Quality Index SQI
%       and extract RR intervals from ECG signal 
%       (single lead ECG signal)
%
%   INPUT:
%       HRVparam : struct of settings for hrv_toolbox analysis
%
%   OUTPUT:
%       rr    :  Vector containing RR interval
%       t     :  Time indices of the rr interval data (seconds)
%       sqijw :  Signal Quality Index
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Function written by Giulia Da Poian
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if nargin < 3
    error('Wrong number of arguments in ConvertRawDataToRRIntervals')
end


if size(ECG_RawData,1)<size(ECG_RawData,2)
    ECG_RawData = ECG_RawData';
end

ECG_RawData = ECG_RawData(:,1); % If more the one leads use only the first

% QRS Dection 1 - jqrs
jqrs_ann = jqrs(ECG_RawData,HRVparams);

% QRS Detection 2 - sqrs
sqrs = run_sqrs(ECG_RawData*2000,HRVparams,0);

% QRS Detection 3 - wqrs
wqrs = wqrsm(ECG_RawData*2000);

% QRS SQI analysis
sqijs = bsqi(jqrs_ann(:),sqrs(:),HRVparams);
sqijw = bsqi(jqrs_ann(:),wqrs(:),HRVparams);

% At this point one can perfrom also Ventricular Fibrillation (VF) 
% and/or Ventricular Tachycardia (VT) analysis


% Translate anotations to rr intervals
rr = diff(jqrs_ann./HRVparams.Fs);
t = jqrs_ann(1:end-1)./HRVparams.Fs;

%%  Export Annotations as ATR files

% Create a Folder for Annotations
WriteAnnotationFolder = [HRVparams.writedata filesep 'Annotation'];
if ~exist(WriteAnnotationFolder, 'dir')
   mkdir(WriteAnnotationFolder)
   fprintf('Creating a new folder: "Annotation", folder is located in %s \n',[pwd filesep WriteAnnotationFolder]);
end
addpath(WriteAnnotationFolder)

% Header File
write_hea([WriteAnnotationFolder filesep subjectID] , HRVparams.Fs, length(ECG_RawData), 'jqrs', 1, 0,'mV')
% ECG QRS and SQI
write_ann([WriteAnnotationFolder filesep subjectID],HRVparams,'jqrs',jqrs_ann);
write_ann([WriteAnnotationFolder filesep subjectID],HRVparams,'sqrs',sqrs);
write_ann([WriteAnnotationFolder filesep subjectID],HRVparams,'wqrs',wqrs);
write_ann([WriteAnnotationFolder filesep subjectID],HRVparams,'sqijs',sqijs);
write_ann([WriteAnnotationFolder filesep subjectID],HRVparams,'sqijw',sqijw);



