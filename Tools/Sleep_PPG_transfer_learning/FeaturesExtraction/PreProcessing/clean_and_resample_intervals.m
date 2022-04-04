function [RR_e_resamp,QRS_A_resamp] = clean_and_resample_intervals(onsets,QRS_A1,Fs)   

% Use same pre-processing used for ECG in the paper ''
% 
% 

IBI = diff(onsets)./Fs;
t_IBI = onsets(2:end)./Fs;

RR_e0 = IBI;     % no ectopic rejection
QRS_A0 = QRS_A1;

% reject ectopic beats
% remove invalid R-R intervals if RR<0.3 or RR>2.0 & 20% outside of 41 points moving mean exclude  

RR_mean_41 = meanfilt1(IBI,41);
RR = [];    % RR interval in seconds
RR_time=[]; % time of RR interval in seconds
v_RR = ones(1,length(IBI))*-1;
n_RR = 0;
for j=1:length(IBI)
    if abs(IBI(j) - RR_mean_41(j)) < 0.2*RR_mean_41(j) &&  IBI(j)>= 0.3 && IBI(j) <= 2.0
        n_RR = n_RR+1;
        RR(n_RR) = IBI(j);
        RR_time(n_RR) = onsets(j)/Fs;
        v_RR(j)=1;
    end
end
valid_RR=find(v_RR==1);
% remove invalid QRS peaks if 50% outside of 41 points moving mean
QRS_A_mean_41=meanfilt1(QRS_A1,41);
n_QRS=0;
QRS=[];
QRS_time=[];
v_QRS=ones(1,length(QRS_A1))*-1;
for k=1:length(QRS_A1)
    if abs(QRS_A1(k)-QRS_A_mean_41(k))<0.5*QRS_A_mean_41(k)
        n_QRS=n_QRS+1;
        QRS(n_QRS) = QRS_A1(k);
        QRS_time(n_QRS) = tIBI(j);
        v_QRS(k)=1;
    end
end
valid_QRS=find(v_QRS==1);
valid_RRQRS=intersect(valid_RR,valid_QRS);
valid_RR=valid_RRQRS;

onsets = onsets(valid_RR);
QRS_A1=QRS_A1(valid_RR);

RR_e1 = IBI(valid_RR);

% resample to 4Hz
onsets = onsets./Fs;
t2=round(onsets(end));
Fs_resamp=4;
t1=(1/Fs_resamp):(1/Fs_resamp):t2; % 4Hz resample

% linear interpolation
QRS_A_resamp = interp1(onsets,QRS_A1,t1);%,'spline');
RR_e_resamp = interp1(onsets,RR_e1,t1);%,'spline');
QRS_A_resamp(find(isnan(QRS_A_resamp)))=nanmean(QRS_A_resamp);
RR_e_resamp(find(isnan(RR_e_resamp)))=nanmean(RR_e_resamp);