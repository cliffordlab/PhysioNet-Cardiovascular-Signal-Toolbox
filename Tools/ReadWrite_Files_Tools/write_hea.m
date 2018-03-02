function write_hea(recordName, fs, datapoints, annotator, gain, offset,unit)

% write_hea(recordName, fs, datapoints, annotator, gain, offset,unit)
%
% ORIGINAL SOURCE AND AUTHORS:     
%       This script written by Qiao Li  
% COPYRIGHT (C) 2016 
% LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

numsig = 1;
%fs = 125;
%datapoints = 75000;
%recordName = '100';
filename = strcat(recordName, '.', annotator);
unit = 'mV';

fileID = fopen(strcat(recordName, '.hea'),'w');
fprintf(fileID,'%s %d %d %d\n', recordName, numsig, fs, datapoints);
for i = 1:numsig
    fprintf(fileID,'%s 16+24 %d/%s 12\n', filename, gain,unit);
end
fprintf(fileID,'#Creator: HRV_toolbox write_hea.m');
fclose(fileID);


end