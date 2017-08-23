function write_hea(recordName, fs, datapoints, annotator, gain, offset,unit)
numsig = 1;
%fs = 125;
%datapoints = 75000;
%recordName = '100';
filename = [recordName '.' annotator];
unit = 'mV';

fileID = fopen([recordName '.hea'],'w');
fprintf(fileID,'%s %d %d %d\n', recordName, numsig, fs, datapoints);
for i = 1:numsig
    fprintf(fileID,'%s 16+24 %d/%s 12\n', filename, gain,unit);
end
fprintf(fileID,'#Creator: HRV_toolbox write_hea.m');
fclose(fileID);


end