fid=fopen('../slpdb/RECORDS','r');
cc=textscan(fid,'%s');
cc=cc{1};
fclose(fid);

for i=1:length(cc)
% read ECG data from SLPDB
[tm,signal,Fs,siginfo]=rdmat(['../slpdb/' cc{i} 'm']);
Fs=Fs(1);

description=squeeze(struct2cell(siginfo));
channel=description(4,:); % or 4, 5, 9 based on rdmat version
ecg_ind=get_index(channel,{'ECG'});
ecg_data=signal(1:end,ecg_ind(1));
ecg_data(find(isnan(ecg_data)))=nanmean(ecg_data);
clear signal

% call sleep staging prediction
class = 4;
predict_label = do_sleep_staging_prediction(ecg_data,Fs,class);
save([cc{i} '_predict_' num2str(class)],'predict_label');
end
    
