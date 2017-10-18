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
%   Written by Giulia Da Poian (giulia.dap@gmail.com) 



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









