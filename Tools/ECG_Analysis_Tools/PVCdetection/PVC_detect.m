function PVCs = PVC_detect(signal,sigName,HRVparams)
%
% PVC_detect(signal,Fs,sigName,OutputFolder,th)
% INPUTS: 
%        signal       : Nx1 raw ech signal (in mV)
%        sigName      : recording name
%        HRVparams    : struct of settings for hrv_toolbox analysis

% OUTPUTS : 
%         PVCs        : locations of detected PVC beats, in samples 
%
%   DEPENDENCIES & LIBRARIES:
%       PhysioNet Cardiovascular Signal Toolbox
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   REFERENCE: 
%   Vest et al. "An Open Source Benchmarked Toolbox for Cardiovascular 
%   Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox%   
%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Qiao Li 
%       Dependent scripts written by various authors 
%       (see functions for details)      
%
%    Modified by Giulia Da Poian on May 2018 to update for compatibility 
%    with the toolbox
%
%
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%
%   Mdified on November 8 2018 by Giulia Da Poian: replace input paramters
%   with toolbox HRVparams struct
%
qrs_times=[];
pvc_outputs=[];



if nargin<3
    error('Wrong number of input parameters: see help for more information')
end

Fs = HRVparams.Fs;          % sampling freqency of record
th = HRVparams.PVC.qrsth;   % threshold for QRS detection (default: 0.1)

LengthSamples =  length(signal);

% split long data to one hour each
t_hours = LengthSamples/(Fs*60*60);

if rem( LengthSamples,(Fs*60*60)) > 0    % if not exactly x hours
    remflag = 1;
else
    remflag = 0;
end

overlap = 10; % Adding 10 seconds to the front of each hour, and throwing away the last 5 and first 5 of each hour for edge effects
overlap_2 = overlap/2;


if t_hours >= 1 
    for j=1:t_hours
        s_start= (j-1) * 200 *60 *60 + 1;
        signal(find(isnan(signal)))=0;

        % edge processing
        if j==1
            edge_signal = ones(overlap_2*Fs,1)*median(signal); % we add half of overlap data before the first hour data
        end
        signal = [edge_signal; signal];

        % call detectpvc()
        [qrs_time, pvc_output] = detectpvc2(signal,Fs,th);

        % edge processing
        valid_time=intersect(find(qrs_time>overlap_2*Fs),find(qrs_time<=length(signal)-overlap_2*Fs));
        qrs_time=qrs_time(valid_time)-overlap_2*Fs;
        pvc_output=pvc_output(valid_time);
        edge_signal=signal(end-overlap*Fs+1:end);

        if j==1
            qrs_time=qrs_time+ (s_start-1) ;  % correct for start of signal to current segment offset
        else
            qrs_time=qrs_time+ (s_start-1) - overlap_2*Fs ;  % correct for start of signal to current segment offset
        end
        qrs_times= ([qrs_times qrs_time]);
        pvc_outputs= ([pvc_outputs pvc_output]);
    end
end

if remflag == 1  % now do remainder of signal
    s_start= floor(t_hours) * Fs *60 *60 + 1;
    signal(find(isnan(signal)))=0;

    % edge processing
    if t_hours<1
        edge_signal=ones(overlap_2*Fs,1)*median(signal); % we add half of overlap data before the first hour data
    end
    signal=[edge_signal; signal];

    % call detectpvc()
    [qrs_time, pvc_output] = detectpvc2(signal,Fs,th);

    % edge processing
    valid_time=find(qrs_time>overlap_2*Fs);
    qrs_time=qrs_time(valid_time)-overlap_2*Fs;
    pvc_output=pvc_output(valid_time);
    
    if t_hours<1
        qrs_time=qrs_time+ (s_start-1) ;  % correct for start of signal to current segment offset
    else
        qrs_time=qrs_time+ (s_start-1) - overlap_2*Fs ;  % correct for start of signal to current segment offset
    end
    
    
    qrs_times= ([qrs_times qrs_time]);
    pvc_outputs= ([pvc_outputs pvc_output]);
end


PVCs = qrs_times(find(pvc_outputs==1));

% write output
annType = repmat('N',length(qrs_times),1);
annType(find(pvc_outputs==1)) = 'V';

% Create a Folder for Annotations
AnnDir = [HRVparams.writedata filesep 'Annotation'];
if ~exist(AnnDir, 'dir')
   mkdir(AnnDir)
   fprintf('Creating a new folder: "Annotation", folder is located in %s \n',[pwd filesep AnnDir]);
end
addpath(AnnDir)
write_ann([AnnDir filesep sigName],HRVparams,'pvc',qrs_times,annType);

