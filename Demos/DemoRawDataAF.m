%	OVERVIEW:
%       This demonstration analyzes a segment of 5-minutes 'raw' data  
%       with known atrial fibrillation to show the operation of the 
%       AF detection algorithm.
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
clear; clc; close all;

HRVparams = InitializeHRVparams('demoAF'); % include the project name
HRVparams.poincare.on = 0; % Poincare analysis off for this demo
HRVparams.DFA.on = 0; % DFA analysis off for this demo

[subjectIDs,filesTBA] = GenerateListOfFilesTBA(HRVparams.ext,HRVparams.readdata,0);
idx = find(strcmp(subjectIDs,'TestAFdata'));
i_patient = idx;

% 1. Load Raw Patient Data

load(filesTBA{i_patient});
% 2. Analyze data using HRV VOSIM toolbox
[results, resFilenameHRV] = Main_VOSIM(signal(:,1),[],'ECGWaveform',HRVparams,subjectIDs(i_patient));


% 3. Compare generated output file with the reference one
        
currentFile = strcat(HRVparams.writedata, filesep, resFilenameHRV, '.csv');
referenceFile = strcat('ReferenceOutput', filesep, 'AFDemo.csv');
testHRV = CompareOutput(currentFile,referenceFile);

% 3. Load QRS annotation saved by Main_VOSIM 
annotName = strcat(HRVparams.writedata, filesep, 'Annotation',filesep,subjectIDs(i_patient));
jqrs_ann = read_ann( annotName{1} , 'jqrs');
wqrs_ann = read_ann( annotName{1} , 'wqrs');

% For demo pourpose recompute bsqi
[sqijw, StartIdxSQIwindows] = bsqi(jqrs_ann,wqrs_ann,HRVparams);

HRVparams.gen_figs = 1;
% Plot detected beats
if HRVparams.gen_figs
    % ECG
    figure;
    plot(time,signal(:,1)); hold on;
    stairs(StartIdxSQIwindows,sqijw);
    for j = 1:length(jqrs_ann)
        line([jqrs_ann(j)./HRVparams.Fs jqrs_ann(j)./HRVparams.Fs],[-1 2],'Color','green');
    end
	legend('ECG','jqrs');
    xlabel('Time (s)'); ylabel('Amplitude (mV)');
    title('ECG')
end


if testHRV 
    fprintf('** DemoRawDataAF: TEST SUCCEEDED ** \n ')
     fprintf('A file named %s.csv \n has been saved in %s \n', ...
    resFilenameHRV, HRVparams.writedata);
else
    fprintf('** DemoRawDataAF: TEST FAILED ** \n')
end


