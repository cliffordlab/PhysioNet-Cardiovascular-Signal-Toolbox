%	OVERVIEW:
%       Demo using 'raw' 5-minutes long ECG signal recoded during  
%       atrial fibrillation (AF) episodes.%        
%   OUTPUT:
%       Detected AF and HRV Metrics 
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Main script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

HRVparams = InitializeHRVparams('rawdatademo'); % include the project name
% Check existence of Input\Output data folders and add to search path
if  isempty(HRVparams.readdata) || ~exist([pwd filesep HRVparams.readdata], 'dir')    
    error('Invalid data INPUT folder');    % If folder name is empty
end
addpath(HRVparams.readdata)


addpath(HRVparams.readdata)
HRVparams.writedata = [HRVparams.writedata filesep 'AF'];
if ~exist(HRVparams.writedata, 'dir')
   mkdir(HRVparams.writedata)
end
addpath(HRVparams.writedata)

[subjectIDs,filesTBA] = GenerateListOfFilesTBA(HRVparams.ext,HRVparams.readdata,0);

idx = find(strcmp(subjectIDs,'TestAFdata'));

i_patient = idx;

%% 1. Process Raw Patient Data

load(filesTBA{i_patient});

HRVparams.Fs = 128;

% QRS Dection 1 - jqrs
[jqrs_ann,sign,en_thres] = jqrs(signal(:,1),HRVparams);

% QRS Detection 2 - sqrs
[sqrs] = run_sqrs(signal(:,1)*2000,HRVparams,0);

% QRS Detection 3 - wqrs
[wqrs]=wqrsm(signal(:,1)*2000);

% QRS SQI
[sqijs,~] = bsqi(jqrs_ann,sqrs,HRVparams);  

[sqijw,RRwindowStartIndices] = bsqi(jqrs_ann,wqrs,HRVparams);

HRVparams.gen_figs = 1;
if HRVparams.gen_figs
    % ECG
    figure;
    plot(time,signal(:,1)); hold on;
    stairs(RRwindowStartIndices,sqijw);
    for j = 1:length(jqrs_ann)
        line([jqrs_ann(j)./HRVparams.Fs jqrs_ann(j)./HRVparams.Fs],[-1 2],'Color','red');
    end
	legend('ECG','jqrs','sqrs','wqrs','SQI JvW','SQI JvS');
    xlabel('Time (s)'); ylabel('Amplitude (mV)');
    title('ECG')
end

%% 2. Export Annotations as ATR files

AnnotationFolder = [HRVparams.writedata filesep 'Annotation' filesep];
if ~exist(AnnotationFolder)
   mkdir(AnnotationFolder);
end
addpath(AnnotationFolder);

% Header File
write_hea([AnnotationFolder subjectIDs{i_patient}], HRVparams.Fs, length(signal), 'jqrs', 1, 0,'mV')

% ECG
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'jqrs',jqrs_ann);
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'sqrs',sqrs);
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'wqrs',wqrs);
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'sqijs',sqijs);
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'sqijw',sqijw);


%% 3. Preprocess RR Data - Using HRV Toolbox
% Remove noise, Remove ectopy, Don't detrend (yet)

rr = diff(jqrs_ann./HRVparams.Fs);
t = jqrs_ann(1:end-1)./HRVparams.Fs;

[NN, tNN, fbeats] = RRIntervalPreprocess(rr,t,[], HRVparams);

%% 4. Calculate Windows
RRwindowStartIndices = CreateWindowRRintervals(tNN, NN, HRVparams);

%% 5. Calculate AF Features

[AFtest, AFwindowsStartIndices,AFfile] = PerformAFdetection(subjectIDs{i_patient},tNN,NN,HRVparams);  
figure(1);
graphannot(AFtest, AFwindowsStartIndices,.25);  

%% 7. Calculate time domain HRV metrics - Using HRV Toolbox
[NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt, SDNN, NNiqr, ...
    RMSSD,pnn50,btsdet,avgsqi,fbeatw, RRwindowStartIndices] = ...
    EvalTimeDomainHRVstats(NN,tNN,[],HRVparams,RRwindowStartIndices);

%% 8. Frequency domain HRV metrics (LF HF TotPow)
%       All Inputs in Seconds

[ulf, vlf, lf, hf, lfhf, ttlpwr, methods, fdflag, window] = ...
   EvalFrequencyDomainHRVstats(NN,tNN, [],HRVparams,RRwindowStartIndices);

%% 9. PRSA
try
    [ac,dc,~] = prsa(NN, tNN, [], RRwindowStartIndices, HRVparams);
catch
    ac = NaN; 
    dc = NaN;
    error_flag(i_patient) = subjectIDs(i_patient);
end

%% 10. SDANN and SDNNi
[SDANN, SDNNI] = ClalcSDANN(RRwindowStartIndices, tNN, NN(:),HRVparams);

%% 12. Export HRV Metrics as CSV File
%Uncomment the following lines for All Results
results = [RRwindowStartIndices(:), ac(:),dc(:),...
    ulf(:),vlf(:),lf(:),hf(:),lfhf(:),ttlpwr(:),fdflag(:),...
    NNmean(:),NNmedian(:),NNmode(:),...
    NNvariance(:),NNskew(:),NNkurt(:),SDNN(:),NNiqr(:),RMSSD(:),pnn50(:),btsdet(:),fbeatw(:)];

col_titles = {'t_win','ac','dc','ulf','vlf','lf','hf','lfhf',...
    'ttlpwr','fdflag','NNmean','NNmedian','NNmode','NNvar','NNskew',...
    'NNkurt','SDNN','NNiqr','RMSSD','pnn50','beatsdetected','corrected_beats'};


% Generates Output - Never comment out
resFilenameHRV = GenerateHRVresultsOutput(subjectIDs{i_patient}, ...
    RRwindowStartIndices,results,col_titles, [],HRVparams, tNN, NN);



%% 11 Compare generated output file with the reference one
        
currentFile = [HRVparams.writedata filesep AFfile '.csv'];
referenceFile = ['ReferenceOutput' filesep 'AFDemo.csv'];
testHRV = CompareOutput(currentFile,referenceFile);

if testHRV 
    fprintf('** DemoRawDataAF: TEST SUCCEEDED ** \n ')
     fprintf('A file named %s.csv \n has been saved in %s \n', ...
    AFfile, HRVparams.writedata);
else
    fprintf('** DemoRawDataAF: TEST FAILED ** \n')
end


