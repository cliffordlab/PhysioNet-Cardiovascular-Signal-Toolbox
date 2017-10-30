%	OVERVIEW:
%       This function demonstrates the function of the synthetic RR interval 
%       generator RRGEN and the calculation of HRV metrics
%   INPUT:
%       No input necessary. RRGEN data is generated within this script
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
%	COPYRIGHT (C) 2017
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
HRVparams = InitializeHRVparams('demo_RRgen'); % include the project name
HRVparams.poincare.on = 0;

% 1. Generate Data

rr = rrgen(HRVparams.demo.length,HRVparams.demo.pe,HRVparams.demo.pn,HRVparams.demo.seed);
t = cumsum(rr);

% 2. Analyze RR Data - Using HRV Toolbox Main function

[results, resFilenameHRV] = Main_VOSIM(rr,t,'RRIntervals',HRVparams,'rrgenData',[]);


% 3 Compare generated output file with the reference one
        
currentFile = [HRVparams.writedata filesep resFilenameHRV '.csv'];
referenceFile = ['ReferenceOutput' filesep 'StandardizedData_HRV_allwindows.csv'];
testHRV = CompareOutput(currentFile,referenceFile);

if testHRV
    fprintf('** DemoStandardizedData: TEST SUCCEEDED ** \n ')
     fprintf('A file named %s.csv \n has been saved in %s \n', ...
    resFilenameHRV, HRVparams.writedata);
else
    fprintf('** DemoStandardizedData: TEST FAILED ** \n')
end
