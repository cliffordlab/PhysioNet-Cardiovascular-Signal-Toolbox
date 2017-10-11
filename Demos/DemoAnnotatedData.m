%	OVERVIEW:
%       Demo using annotated data 
%   OUTPUT:
%       HRV Metrics
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
HRVparams = InitializeHRVparams('demo');   % Initialize settings for demo
% Check existence of Input\Output data folders and add to search path

if  isempty(HRVparams.readdata) || ~exist([pwd filesep HRVparams.readdata], 'dir')    
    error('Invalid data INPUT folder');    % If folder name is empty
end
addpath(HRVparams.readdata)


addpath(HRVparams.readdata)
HRVparams.writedata = [HRVparams.writedata filesep 'Annotated'];
if ~exist(HRVparams.writedata, 'dir')
   mkdir(HRVparams.writedata)
end
addpath(HRVparams.writedata)


% Check for a list of files to be analyzed in current directory
% in .mat format

[subjectIDs,filesTBA] = GenerateListOfFilesTBA(HRVparams.ext, HRVparams.readdata,[]);

% Prepare for parallel loop by eliminating variables
clear nummatchingfiles x i filename flag match
numsub = length(subjectIDs);

% NOTE: This loop can be run in parallel by changing the loop to a parfor
% loop.

for i_patient = 1:numsub   
    try
        
        %% 1. Import Patient Data
        % HRVparams = initialize_HRVparams('demo'); % enable this when using parfor loops
        RRwindowStartIndices = [];
        tNN = [];
        NN = [];
        [samples, t, rr, ~, ~, annotations, HRVparams.Fs, ~, ~] = read_qrs(filesTBA{i_patient},HRVparams.datatype);
        
        %% 2. Perform HRV analysis on the RR intervals
        HRVparams.MSE.on = 0; % No MSE analysis for this demo
        [results, resFilenameHRV] = Main_VOSIM(rr,t,'RRIntervals',HRVparams,subjectIDs(i_patient),annotations);

              
        %% 3. Compare generated output file with the reference one
        
        currentFile = [HRVparams.writedata filesep resFilenameHRV '.csv'];
        referenceFile = ['ReferenceOutput' filesep 'Annotated_HRV_allwindows.csv'];
        testHRV = CompareOutput(currentFile,referenceFile);
        
        if testHRV
            fprintf('** DemoAnnotatedData: TEST SUCCEEDED ** \n ')
        else
            fprintf('** DemoAnnotatedData: TEST FAILED ** \n')
        end
        
    catch
        if isnumeric(subjectIDs(i_patient))
            current_filename = ['error_' num2str(subjectIDs(i_patient))];
        else
            current_filename = ['error_' subjectIDs(i_patient)];
        end
            
        results = NaN;
        col_titles = {'NaN'};
        fprintf('Error on subject %s \n', char(subjectIDs(i_patient)));

        fprintf('** DemoAnnotatedData: TEST FAILED ** \n')
    end
    
end



