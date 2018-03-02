function [samples, time, rr_sec, rr_samp, HR, annotations, Fs, tstart, date] = read_qrs(record,datatype)
% [samples,time,rr_sec,rr_samp,HR,annotations,Fs] = read_qrs('1003.qrs')
%	OVERVIEW:   Reads the GE QRSDK format QRS detection 
%               file and converts to time-in-samples and annotations
%   INPUT:      record - QRSDK RR interval file name in string format
%               datatype - type of data being read into function
%                   'MARS' - default!
%                   THE FOLLOWING OPTIONS ARE NOT OPERATIONAL YET!
%                   'CSV'
%                   'WFDB'
%   OUTPUT:     samples - location of QRS annotations in samples from the 
%                   beginning of the data
%               time - location of QRS ann in time (seconds)
%               rr_sec - rr interval values in seconds
%               rr_samp - rr interval values in samples
%               HR - beat-by-beat Heart Rate
%               annotations - QRS annotation labels
%                   ? N = Normal
%                   ? B = Bundle Branch Block
%                   ? A = Aberrant
%                   ? V = Ventricular
%                   ? F = Ventricular Fusion
%                   ? S = Supraventricular
%                   ? E = Ventricular Esc.
%                   ? J = Junctional
%                   ? P = Ventricular Paced
%                   ? I = Idioventricular
%                   ? 00 = no signal
%                   ? 01 = noise on channel 1
%                   ? 02 = noise on channel 2
%                   ? 03 = noise on both channels (unreadable)
%               Fs - sampling frequency
%
%   ORIGINAL SOURCE AND AUTHORS:     
%       Qiao Li 
%       
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%
if nargin < 1
    error('No input supplied')
end
if nargin < 2
    datatype = 'MARS';
end

fid=fopen(record,'r');

% Adding the following for working on LINUX:

% move the pointer of file to offset 60 for reading starting time
% 8 characters + 2 NULL terminators
status=fseek(fid, 60, -1);
tstart=fread(fid,8,'*char'); % ECG Start Time HH:MM:SS
% move the pointer of file to offset 70 for reading starting date
status=fseek(fid, 70, -1);
date=fread(fid,9,'*char');%ECG Start Date DD-MMM-YY
% move the pointer of file to offset 90 for reading Fs
status=fseek(fid, 90, -1);
Fs=fread(fid,1,'short');
% timeincrement = 1/Fs;
% move the pointer of file to the offset 512 for annotations
status=fseek(fid,512,-1);
[ann_pair count]=fread(fid,2,'uint8');
ann_number=0;
time_accumulation=0;
while (count==2)
    time=ann_pair(1);
    ann=ann_pair(2);
    while (ann_pair(2)==0 & ann_pair(1)==255) % if more than 255 samples between beats, read additional blocks
        [ann_pair count]=fread(fid,2,'uint8');
        time=time+ann_pair(1);
        ann=ann_pair(2);
    end
    time_accumulation=time_accumulation+time;
    ann_number=ann_number+1;
    samples(ann_number)=time_accumulation;
    annotations{ann_number,1}=char(ann);
    rr_samp(ann_number)=time;
    [ann_pair count]=fread(fid,2,'uint8');
end
rr_samp(1) = []; % delete first point
rr_samp(length(rr_samp)) = [];% delete last point
samples(1) = [];
samples(length(samples)) = [];% delete last point
annotations(1) = [];
annotations(length(annotations)) = [];% delete last point

HR = 60.0 ./ (rr_samp./Fs);
time = samples./Fs;
rr_sec = rr_samp./Fs; 
% RRTachogram = [time',rr_sec,ann];
fclose(fid);
end
    


