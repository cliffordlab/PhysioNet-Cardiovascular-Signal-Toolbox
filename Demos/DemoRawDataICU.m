%	OVERVIEW:
%       This demonstration analyzes a segment of data collected in the 
%       intensive care unit (ICU) which contains ECG, ABP, and PPG signals 
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

clear; clc; close all;

HRVparams = InitializeHRVparams('demoICU'); % include the project name
HRVparams.poincare.on = 0; % Pinocare analysis off for this demo
HRVparams.DFA.on = 0; % DFA analysis off for this demo


[subjectIDs,filesTBA] = GenerateListOfFilesTBA(HRVparams.ext,HRVparams.readdata,0);

idx = find(strcmp(subjectIDs,'TestICUdata'));

i_patient = idx;

% 1. Load Raw Patient Data (ECG Waveform)
load(filesTBA{i_patient});

% 2. Analyze data using HRV VOSIM toolbox
[~, resFilenameHRV] = Main_VOSIM(signal(:,1),[],'ECGWaveform',HRVparams,subjectIDs(i_patient),[],signal(:,5),'PPG',signal(:,3),'ABP');

% 3. Load annotations ans SQI for plot
AnnFile = strcat(HRVparams.writedata, filesep, 'Annotation', filesep, subjectIDs{i_patient});
jqrs_ann = read_ann(AnnFile,'jqrs');

ppg_ann = rdann(AnnFile,'ppg');
qppg(signal(:,5),HRVparams.Fs);
[~,ppgsqi,~] = read_ann(AnnFile,'sqippg');

abpann = read_ann(AnnFile,'abpm');
features =  abpfeature(signal(:,3), abpann, HRVparams.Fs);
[BeatQ, goodbeats] = jSQI(features, abpann, signal(:,3));

% 5. Plotting
HRVparams.gen_figs = 1;
if HRVparams.gen_figs
    % ECG
    figure;
    for j = 1:length(jqrs_ann)
        line([jqrs_ann(j)./HRVparams.Fs jqrs_ann(j)./HRVparams.Fs],[-1 2],'Color','red');
        %plot(jqrs_ann./s.Fs,.8.*ones(length(jqrs_ann),1),'o'); hold on;
    end
    clear j
    hold on;
	plot(time,signal(:,1)); hold on;
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

% 6. Pulse Transit Time

ptt = pulsetransit(jqrs_ann, abpann);

% 7. Plot BP vs PTT
syst = features(1:length(ptt),2);

if HRVparams.gen_figs
    figure;
    plot(syst,ptt(:,3)./HRVparams.Fs,'o');
    xlabel('BP (mmHg)'); ylabel('PTT (s)');
    title('Pulse Transit Time - BP vs PTT (ABP - QRS)')
	
end


% 8. Compare generated output file with the reference one
        
currentFile = strcat(HRVparams.writedata, filesep, resFilenameHRV, '.csv');
referenceFile = strcat('ReferenceOutput', filesep, 'ICU_HRV_allwindows.csv');
testHRV = CompareOutput(currentFile,referenceFile);

if testHRV
    fprintf('** DemoRawDataICU: TEST SUCCEEDED ** \n ')
     fprintf('A file named %s.csv \n has been saved in %s \n', ...
    resFilenameHRV, HRVparams.writedata);
else
    fprintf('** DemoRawDataICU: TEST FAILED ** \n')
end