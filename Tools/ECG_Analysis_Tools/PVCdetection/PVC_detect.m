function PVC_detect(path,record,FS,wrann_text,th)
% input: path - path of record
%        record - WFDB record name
%        FS - sampling freqency of record, may be overwritten by WFDB 
%           functions
%        wrann_text - write annotation by a text file for long data 
%           (if wrann failed)
%        th - threshold for QRS detection (default: 0.1), if too many
%           missing beats, decrease th; if too many extra beats, 
%           increase th
% version: 7v
%
% Emailed from Qiao Li to Adriana Vest on May 30th 2017
% Modified July 21, 2017 to update for compatibility with the HRV toolbox
% Modifications include all things commented with ANV
%
%

qrs_times=[];
pvc_outputs=[];

if nargin<5
    th=0.1;
end

if nargin<4
    wrann_text=0;
end

if nargin<3
    FS=360;
end

% Create parallel pool 
%poolobj = parpool;   % ANV commented out

display('Read data ...');

[siginfo,Fs] = wfdbdesc([path record]);

% ANV: added the following for records that do not have a header file:
if isempty(siginfo)
    [tm, signal] = rdsamp([path record], [1]); 
    siginfo(1).LengthSamples = length();
end

if isempty(Fs)
    Fs=FS;
end
Fs=Fs(1);

% split long data to one hour each
t_hours=siginfo(1).LengthSamples/(Fs*60*60);
t_days= floor(t_hours/24);
if rem( siginfo(1).LengthSamples,(Fs*60*60)) > 0    % if not exactly x hours
    remflag=1;
else
    remflag=0;
end

overlap=10; % Adding 10 seconds to the front of each hour, and throwing away the last 5 and first 5 of each hour for edge effects
overlap_2=overlap/2;

X = sprintf('Analyzing %d hours of data',t_hours);
disp(X)
signals=[];
s_start=1;
if t_hours >= 1 % ANV added statement to enable processing on records < 1 hour
    parfor j=1:t_hours
    X = sprintf('Analyzing hour number %d',j);
    disp(X)
    s_start= (j-1) * 200 *60 *60 + 1;
    s_end= j * 200 *60 *60;
    % removed Fs from rdsamp return to avoid Fs reading by mistake
    % read the first [1] channel of data
    [tm, signal]=rdsamp([path record], [1], s_end, s_start); 
%     if isempty(Fs)
%         Fs=FS;
%     end

    %signal=signal(:,1);
    signal(find(isnan(signal)))=0;

    % edge processing
    if j==1
        edge_signal=ones(overlap_2*Fs,1)*median(signal); % we add half of overlap data before the first hour data
    end
    signal=[edge_signal; signal];

    % call detectpvc()
    [qrs_time pvc_output] = detectpvc2(signal,Fs,th);

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

if remflag==1  % now do remainder of signal
    s_start= floor(t_hours) * Fs *60 *60 + 1;
    s_end= siginfo(1).LengthSamples-1;
    % removed Fs from rdsamp return to avoid Fs reading by mistake
    [tm,signal] = rdsamp([path '/' record], [1], s_end, s_start); % ANV added '/'

    %signal=signal(:,1);
    signal(find(isnan(signal)))=0;

    % edge processing
    if t_hours<1
        edge_signal=ones(overlap_2*Fs,1)*median(signal); % we add half of overlap data before the first hour data
    end
    signal=[edge_signal; signal];

    % call detectpvc()
    [qrs_time pvc_output] = detectpvc2(signal,Fs,th);

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

display('Writing data ...');
% write output
annType=[];
annType(find(pvc_outputs==0))='N';
annType(find(pvc_outputs==1))='V';
num=zeros(length(qrs_times),1);

if wrann_text==0
    wrann_path([path '/'],record,'pvc_d',qrs_times,annType,'0',0,num); % ANV added '/'
else
    wrann_by_text(Fs,record,'pvc_d',qrs_times,annType,'0',0,num);
end

display('Finished.');

% Delete the parallel pool
%delete(poolobj);
end