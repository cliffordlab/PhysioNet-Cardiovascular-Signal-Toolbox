%	OVERVIEW:
%       This is demo for VOSIM HRV toolbox using RR intervals with
%       annotations. Provided data are a subset from the MIT Physionet 
%       NSR dataset, which contains long-term ECG recordings of subjects 
%       in normal sinus rhythm.
%       It shows how to automaticly import multiple files from a folder, 
%       perfrom the HRV analysis on each of them and then store the results 
%       in .csv format.  
%       It uses the default parameters in the configuration file using 'demo'
%       option : InitializeHRVparams('demo_NSR').
%
%   OUTPUT:
%       HRV Metrics exported to .cvs files
%
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

 % Initialize settings for demo
HRVparams = InitializeHRVparams('demo_NSR');  
HRVparams.MSE.on = 0; % No MSE analysis for this demo
HRVparams.DFA.on = 0; % No DFA analysis for this demo
HRVparams.output.separate = 0;   % For this demo write all the results in one file 

% Check for a list of files to be analyzed in current directory in .dat format
[subjectIDs,filesTBA] = GenerateListOfFilesTBA(HRVparams.ext, HRVparams.readdata,[]);

% Prepare for parallel loop by eliminating variables
clear nummatchingfiles x i filename flag match
numsub = length(subjectIDs);
notAnalyzed = 0;

% NOTE: This loop can be run in parallel by changing the loop to a parfor
% loop.
for i_patient = 1:numsub   
    
    thisPatient = strcat(HRVparams.readdata, filesep, subjectIDs(i_patient));
    try
        
        % 1. Import Patient Data
        RRwindowStartIndices = [];
        tNN = [];
        NN = [];
        [samples,annotations] = rdann(thisPatient{1},HRVparams.ext);
        rr = diff(samples)./HRVparams.Fs; 
        t = cumsum(rr);
        
        % Demo keeps only the first 2h 
        
        rr = rr(t<60*60*2);
        t = t(t<60*60*2);
               
        % 2. Perform HRV analysis on the RR intervals
        [results, resFilenameHRV] = Main_VOSIM(rr,t,'RRIntervals',HRVparams,subjectIDs(i_patient),annotations);
        currentFile = [HRVparams.writedata filesep resFilenameHRV '.csv'];

        
    catch
       
        results = NaN;
        col_titles = {'NaN'};
        currentFile = '';
        notAnalyzed = 1;
        fprintf('Error on subject %s \n', char(subjectIDs(i_patient)));    

    end
    
end


% 3. Compare generated output file with the reference one

referenceFile = ['ReferenceOutput' filesep 'NSR_HRV_allwindows_allpatients.csv'];
testHRV = CompareOutput(currentFile,referenceFile);

if testHRV
    fprintf('** DemoAnnotatedData: TEST SUCCEEDED ** \n ')
elseif notAnalyzed == 0
    fprintf('** DemoAnnotatedData: TEST FAILED ** \n')
    fprintf('Error: generated output does not match reference \n')
elseif notAnalyzed == 1
    fprintf('** DemoAnnotatedData: TEST FAILED ** \n')
    fprintf('Error: analysis not performed \n');    
end
