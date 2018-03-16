function RemovedWindowsStats(Windows, AFWindows, HRVparams, sub_id)

%   RemovedWindowsStats(Windows, AFWindows, HRVparams, sub_id)
%	OVERVIEW:
%       Returns a file with the total %ages of windows removed, the %ages 
%       removed due  to AF and  %ages low quality signal
%
%   INPUT:
%       Windows    : vector containing windows indexes corresponding 
%                    to RR interval segmentation for HRV analysis    
%       AFWindows  : vector containing windows indexes corresponding 
%                    to RR interval segmentation for AF analysiss 
%       HRVparams  : struct containing the parameters for HRV analysis
%       subjectID  : string with subject name
%
%   DEPENDENCIES & LIBRARIES:
%       PhysioNet Cardiovascular Signal Toolbox
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   REFERENCE: 
%   Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%   Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Giulia Da Poian (giulia.dap@gmail.com) on Sep 6, 2017.
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



TotWind = length(Windows);

WindNotAnalyzed = numel(find(isnan(Windows)));

AFWind = numel(Windows(AFWindows));
LowQualityWind = WindNotAnalyzed - AFWind;

percentAF = AFWind/TotWind *100;
percentLowQuality = LowQualityWind/TotWind *100;
percentNotAnalyzed = WindNotAnalyzed/TotWind*100;


filename = strcat(HRVparams.writedata, filesep,'RemovedWindowsStatsSummary_', HRVparams.time,'.csv');
try
     T = readtable(filename);
catch
    T = []; % If file does not exist yet
end


variables_names = [{'patID' 'TotWind' 'NotAnalyzed' 'PercentNotAnalyzed' 'PercentAFWind' 'PercentLowQualityWind'}]; % Add colum with patient ID
patid_array  = string({sub_id});
variables_vals = [patid_array TotWind WindNotAnalyzed percentNotAnalyzed percentAF percentLowQuality]; % Add colum with patient ID

% Write AF results to a table 
T =  [T ; array2table(variables_vals,'VariableNames',variables_names)];
% Use writetable to geberate csv file with the results   
writetable(T,filename);









