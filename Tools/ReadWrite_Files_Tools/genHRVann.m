function genHRVann(record, HRVparams, annName, samples, metrictype, metric)
%
%   genHRVann(record, HRVparams, annName, samples, metrictype, metric)
%
%   OVERVIEW:   
%       This function generates annotation files to be used with mxm.c
%       using the function write_ann.m
%   INPUT:      
%       record     : a string with the name of the subject/record
%       HRVparams  : struct with settings
%       annName    : a string with the desired file extension
%       samples    : the time data of the annotations to be written (in
%                    samples)
%       metrictype : this details the type of measurement to be written to
%                    an annotion file according to the following standard
%                    0 = SDNN      5 = TTLPWR
%                    1 = LF/HF     6 = RMSSD
%                    2 = LF        7 = PNN50
%                    3 = HF     	8 = Mean
%                    4 = VLF    	
%       metric     : the annotation to be saved. 
%
%   OUTPUT:
%       A .annName file with annotation timing and metrics.
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       This script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

type = repmat('=',length(metric),1);  % this type needs to be on every annotation
subtype = repmat(metrictype,length(metric),1);
% this details the type of measurement
                % 0 = SDNN      5 = 
                % 1 = LF/HF     6 = 
                % 2 = RMSSD     7 = 
                % 3 =           8 = 
                % 4 =           9 = 
chan = zeros(length(metric),1);       % Default Value
num =  zeros(length(metric),1);       % Default Value

comments = cell(length(metric));

for k = 1:length(metric)
    try
        comments{k} = char(string(metric(k)));  % Actual HRV Metric to be compared in string format
    catch
        comments{k} = [];
    end
end

write_ann(record,HRVparams,annName,samples,type,subtype,chan,num,comments);

% Convert annotation locations from sample number to WFDB compatible time
%[timeStamp,~] = wfdbtime(record,samples);

% % Write annotations to .txt file in WFDB compatible format
% filename = [record '_temp.txt'];
% for i = 1:length(samples)
%     fileID = fopen(filename,'a');
%     %fprintf(fileID, '\t%s %7d\t%c\t%5d%5d%5d\r\n',time_formatted,samples(i),ann(i),subType(i),chan(i),num(i));
%     fprintf(fileID,'%12s%9d%6s%5d%5d%5d%12s\0\n', timeStamp{i},samples(i),type(i),subtype(i),chan(i),num(i),comments{i});
%     fclose(fileID);
% end
% clear i
% 
% %% 3. Convert .txt file to a WFDB formatted binary .ann file
% statement = ['!wrann -r ' record ' -a ann < ' record '_temp.txt'];
% eval(statement);
% 
% statement = ['!rdann -r ' record ' -a ann > 1test.txt'];
% eval(statement);

end