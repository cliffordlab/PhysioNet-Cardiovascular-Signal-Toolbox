%	OVERVIEW:
%       Demo using ICU rawdata 
%   INPUT:
%       
%   OUTPUT:
%       HRV Metrics
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

HRVparams = InitializeHRVparams('rawdatademo'); % include the project name

HRVparams.writedata = [HRVparams.writedata filesep 'ICU'];
if ~exist([pwd filesep HRVparams.writedata], 'dir')
   mkdir(HRVparams.writedata)
end
addpath(HRVparams.writedata)


[subjectIDs,filesTBA] = GenerateListOfFilesTBA(HRVparams.ext,HRVparams.readdata,0);

idx = find(strcmp(subjectIDs,'TestICUdata'));

i_patient = idx;

%% 1. Process Raw Patient Data

load(filesTBA{i_patient});

% QRS Dection 1 - jqrs
[jqrs_ann,sign,en_thres] = jqrs(signal(:,1),HRVparams);

% QRS Detection 2 - sqrs
[sqrs] = run_sqrs(signal(:,1)*2000,HRVparams,0);

% QRS Detection 3 - wqrs
[wqrs]=wqrsm(signal(:,1)*2000);

% QRS SQI
[sqijs,windows_all] = bsqi(jqrs_ann(:),sqrs(:),HRVparams);

[sqijw,windows_all] = bsqi(jqrs_ann(:),wqrs(:),HRVparams);

% PPG Detection - qppg
[ppg_ann] = qppg(signal(:,5),HRVparams.Fs);

% PPG SQI 
[ppgsqi,~,~,~] = PPG_SQI_buf(signal(:,5),ppg_ann);

% ABP
abpann = run_wabp(signal(:,3));

% ABP SQI
features =  abpfeature(signal(:,3), abpann, HRVparams.Fs);
[BeatQ, goodbeats] = jSQI(features, abpann, signal(:,3));

HRVparams.gen_figs = 1;
if HRVparams.gen_figs
    % ECG
    figure;
    %stairs(windows_all,sqijw); hold on;
    stairs(windows_all,sqijs);
    for j = 1:length(jqrs_ann)
        line([jqrs_ann(j)./HRVparams.Fs jqrs_ann(j)./HRVparams.Fs],[-1 2],'Color','red');
        %plot(jqrs_ann./s.Fs,.8.*ones(length(jqrs_ann),1),'o'); hold on;
    end
    clear j
    hold on;
	plot(time,signal(:,1)); hold on;
    %plot(sqrs./s.Fs,.79.*ones(length(sqrs),1),'o'); hold on;
    %plot(wqrs./s.Fs,.81.*ones(length(wqrs),1),'o'); hold on;
    %plot(windows_all,sqijw,'.','markersize',10); hold on;
    %plot(windows_all, sqijs,'.','markersize',10); hold on;
	%legend('ECG','jqrs','sqrs','wqrs','SQI JvW','SQI JvS');
    xlabel('Time (s)'); ylabel('Amplitude (mV)');
    title('ECG')
    
    % ABP
    figure;
    plot(time,signal(:,3)); hold on;
    for j = 1:length(abpann)
        line([abpann(j)./HRVparams.Fs abpann(j)./HRVparams.Fs],[0 200],'Color','red');
        %plot(jqrs_ann./s.Fs,.8.*ones(length(jqrs_ann),1),'o'); hold on;
    end
    clear j
    %plot(abpann./s.Fs,100.*ones(length(abpann)),'o');
    %graphannot(BeatQ(:,1), abpann./s.Fs,110,s);
    stairs(abpann./HRVparams.Fs,100*double(~BeatQ(:,1)));
    xlabel('Time (s)'); ylabel('Amplitude (mmHg)');
    title('ABP')
    
    % PLETH
    figure;
    plot(time,signal(:,5)); hold on;
	for j = 1:length(ppg_ann)-1
        line([ppg_ann(j)./HRVparams.Fs ppg_ann(j)./HRVparams.Fs],[-1 1],'Color','red');
        %plot(jqrs_ann./s.Fs,.8.*ones(length(jqrs_ann),1),'o'); hold on;
        if (strcmp(ppgsqi(j),'A') || strcmp(ppgsqi(j),'E'))
            numsqi(j) = 1;
        else
            numsqi(j) = 0;
        end
    end
    numsqi(j+1) = 0;
    clear j
	%plot(ppg_ann./s.Fs,ones(length(ppg_ann)),'o');
    %graphannot(ppgsqi,ppg_ann./s.Fs,1.1,s);
    stairs(ppg_ann./HRVparams.Fs,numsqi);
    xlabel('Time (s)'); ylabel('Amplitude (mV)');
    title('PPG')
    
    % RESP
    figure;
    plot(time,signal(:,6));
    title('Respiration')

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
% PPG
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'ppg',ppg_ann);
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'sqippg',ppg_ann,char(ppgsqi));
% ABP
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'abpm',abpann);
write_ann([AnnotationFolder subjectIDs{i_patient}],HRVparams,'sqiabp',BeatQ(:,1));

%% 3. Pulse Transit Time

ptt = pulsetransit(jqrs_ann, abpann);

%% Plot BP vs PTT
syst = features(:,2);

if HRVparams.gen_figs
    figure;
    plot(syst,ptt(:,3)./HRVparams.Fs,'o');
    xlabel('BP (mmHg)'); ylabel('PTT (s)');
    title('Pulse Transit Time - BP vs PTT (ABP - QRS)')
	
end


%% 4. Preprocess RR Data - Using HRV Toolbox

rr = diff(jqrs_ann./HRVparams.Fs);
t = jqrs_ann(1:end-1)./HRVparams.Fs;

% Remove noise, Remove ectopy, Don't detrend (yet)
[NN, tNN, ~] = RRIntervalPreprocess(rr,t,[], [], HRVparams);

%% 5. Calculate Windows
windows_all = CreateWindowRRintervals(tNN, NN, HRVparams);

%% 6. Calculate AF Features

afresults = PerformAFdetection(subjectIDs{i_patient},tNN,NN,HRVparams);

%% 7. Calculate time domain HRV metrics - Using HRV Toolbox
fbeats = zeros(length(NN),1);
[NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt, SDNN, NNiqr, ...
    RMSSD,pnn50,btsdet,avgsqi,fbeatw, windows_all] = ...
    EvalTimeDomainHRVstats(NN,tNN,[],HRVparams,windows_all,fbeats);

%% 8. Frequency domain HRV metrics (LF HF TotPow)
%       All Inputs in Seconds

[ulf, vlf, lf, hf, lfhf, ttlpwr, methods, fdflag, window] = ...
   EvalFrequencyDomainHRVstats(NN,tNN, [],HRVparams,windows_all);

%% 9. PRSA
try
    [ac,dc,~] = prsa(NN, tNN, [], windows_all, HRVparams);
catch
    ac = NaN; 
    dc = NaN;
end

%% 10. SDANN and SDNNi
[SDANN, SDNNI] = ClalcSDANN(windows_all, tNN, NN(:),HRVparams); 

%% 11. Export HRV Metrics as CSV File
%Uncomment the following lines for All Results
results = [windows_all(:), ac(:),dc(:),...
    ulf(:),vlf(:),lf(:),hf(:),lfhf(:),ttlpwr(:),fdflag(:),...
    NNmean(:),NNmedian(:),NNmode(:),...
    NNvariance(:),NNskew(:),NNkurt(:),SDNN(:),NNiqr(:),RMSSD(:),pnn50(:),btsdet(:),fbeatw(:)];

col_titles = {'t_win','ac','dc','ulf','vlf','lf','hf','lfhf',...
    'ttlpwr','fdflag','NNmean','NNmedian','NNmode','NNvar','NNskew',...
    'NNkurt','SDNN','NNiqr','RMSSD','pnn50','beatsdetected','corrected_beats'};

% Uncomment the following lines for just two HRV metrics
%results = [NNmean(:), NNmedian(:)];
%col_titles = {'NN Mean','NNmedian'};

% Generates Output - Never comment out
resFilename = GenerateHRVresultsOutput(subjectIDs(i_patient),windows_all,results,col_titles, [],HRVparams, tNN, NN);

fprintf('A file named %s.%s \n has been saved in %s \n', ...
    resFilename,HRVparams.output.format, HRVparams.writedata);


fprintf('ICU demo completed with success\n')