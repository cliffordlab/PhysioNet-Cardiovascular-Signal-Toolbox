% demo_ECG_sleep_staging_slpdb
% Note: Raw slpdb data should be converted to Matlab format by wfdb2mat function before analysis.

fid=fopen('RECORDS');
cc=textscan(fid,'%s');
cc=cc{1};
fclose(fid);

class=1;
% class = 1, four types classification
% class = 4, two types classification
%   sleep_stage:  = 1: Wake, 
%                   2: REM sleep, 
%                   3: NREM light sleep, 
%                   4: NREM deep sleep, when class = 1
%   sleep_stage:  = 1: NREM sleep, 
%                   2: Wake+REM, when class = 4

% par
for i=1:length(cc)
    do_ECG_sleep_staging_slpdb_sub(cc{i},class)
end

function do_ECG_sleep_staging_slpdb_sub(record,class)

record
[tm,signal,Fs,siginfo]=rdmat([record 'm']);
ECGdata=signal(:,1);
fs=250;
% class=1;

sleep_stage = ECG_sleep_staging_prediction(ECGdata,fs,class);

% read stage type from annotation
[ann,type,subtype,chan,num,comments]=read_ann(record,'st');
ann_epoch_len=30; % annotation epoch is 30s
seg_length=300; % segment length 300s (5 min)
stage=[];
stage(1:round(ann(end)/Fs/ann_epoch_len)+1)=NaN;

for j=1:length(ann)
    time_j=round(ann(j)/Fs/ann_epoch_len)+1;
    switch (comments{j,1}(1))
        case 'W'
            stage(time_j)=1;
        case 'R'
            stage(time_j)=2;
        case '1'
            stage(time_j)=3;
        case '2'
            stage(time_j)=3;
        case '3'
            stage(time_j)=4;
        case '4'
            stage(time_j)=4;
    end
end
if class==4
    stage(find(stage==1))=2;
    stage(find(stage>2))=1;
end
save(['ECG_sleep_staging_result_class_' num2str(class) '_' record],'sleep_stage','stage');
end