function [NN_resamp,RIAV_resamp,sqi_resamp,t1] = clean_resample_timeseries(onsets,RIAV,sqi,Fs,Fs_resamp)

% Use time-series of peaks and onsets to derive a new time-series which 
% represents information related to respiration, in particular respiratory 
% induced  amplitude  variation
%
% Inputs :
%     onsets : detected PPG onsets
%       RIAV : original respiratory-induced amplitude variation time series
%        sqi : signal quality index of PPG
%         Fs : PPG sampling frequency
%  Fs_resamp : resampling frequency (suggested 4Hz)
%
% Outputs :
%     NN_resamp : clean time series of IBI
%   RIAV_resamp : clean and resampled respiratory-induced amplitude
%                variation time series
%
% Written by Giulia Da Poian, 07-Dec-2018
%


N_medianFilt = 41;

IBI = diff(onsets)./Fs; % Intra beat intervals 
RIAV = RIAV(2:end);       % respiratory amplitude 
t_IBI = onsets(2:end);
sqi = sqi(2:end);

% reject ectopic beats
% remove invalid intervals if RR<0.3 or RR>2.0 & 20% outside of 41 points moving mean exclude  

RR_filt = meanfilt1(IBI,N_medianFilt);

NN = []; % NN intervals in seconds
tNN = []; % time of NN interval in seconds
v_RR = ones(1,length(IBI))*-1;
n_RR=0;

for j=1:length(IBI)
    if abs(IBI(j)-RR_filt(j)) < 0.2*RR_filt(j) && IBI(j)>= 0.3 && IBI(j) <= 2.0
        n_RR=n_RR+1;
        NN(n_RR) = IBI(j);
        tNN(n_RR) = t_IBI(j)/Fs;
        v_RR(j)=1;
    end
end
valid_RR = find(v_RR==1);

% remove invalid peaks if 50% outside of 41 points moving mean
RIAV_filt = meanfilt1(RIAV,N_medianFilt);
n_peaks = 0;
ppgPeak = [];
t_ppgPeak = [];
v_ppg=ones(1,length(RIAV))*-1;
for k=1:length(RIAV)
    if abs(RIAV(k)-RIAV_filt(k))<0.5*RIAV_filt(k)
        n_peaks = n_peaks+1;
        ppgPeak(n_peaks) = RIAV(k);
        t_ppgPeak(n_peaks) = t_IBI(k)/Fs;
        v_ppg(k)=1;
    end
end
valid_ppgPeak = find(v_ppg==1);
valid_RR = intersect(valid_RR,valid_ppgPeak);

RIAV = RIAV(valid_RR);
NN = IBI(valid_RR);
t_NN = t_IBI(valid_RR);

t_NN = t_NN./Fs;
t2 = round(t_NN(end));
t1 = (1/Fs_resamp):(1/Fs_resamp):t2; 

% for sqi
t_sqi = t_IBI./Fs;
t2 = round(t_sqi(end));
tsqi1 = (1/Fs_resamp):(1/Fs_resamp):t2; 

% linear interpolation
RIAV_resamp = interp1(t_NN,RIAV,t1);
NN_resamp = interp1(t_NN,NN,t1);
sqi_resamp = interp1(t_sqi,sqi,tsqi1);
RIAV_resamp(isnan(RIAV_resamp)) = nanmean(RIAV_resamp);
NN_resamp(isnan(NN_resamp)) = nanmean(NN_resamp);
sqi_resamp(isnan(sqi_resamp)) = nanmean(sqi_resamp);
sqi_resamp=sqi_resamp(1:length(NN_resamp));