function imdb = Extract_PPG_feautures(PPG,HRVparams,SigName)

% function imdb = PPG_feat_ext(PPG,HRVparams)
% This function is a modified version of CRC_features_calculation used in
% CV_SleepStaging designed to dela with PPG signals insted ECG
% 
% Inputs:
%        PPG : vector containing the PPG signals
%  HRVparams : structur with 
%    SigName : name use to save derived features
%
% Output:
%       imdb : structur containing the extracted feature to be used for
%              classifictaion
%
% This function is a part of the PPG_sleepstagin project
%         https://github.com/cliffordlab/PPG-sleepstaging
%
% A list of dependencies is given in the readme file
%
% Author : Giulia Da Poian (giulia.dap@gmail.com), 07-Dec-2018
% Part of this code was originally included in the CRC_features_ectract
% in the CV_SleepStaging code, and has been modified to work with PPG
% signals
%


imdb=[];
imdb.images.data=[];
imdb.meta.features=[];  

SavePath = HRVparams.writedata;
Fs_resamp = 4; 
NFFT = 2000;

% 512 points for 2 Hz (Fs_resamp/2)
very_low_band = 1:3;  % <0.01
low_band = 4:25;      % 0.01-0.1
high_band = 26:103;   % 0.1-0.4

win_inc = HRVparams.increment; % sliding window length is 30s
win_len = HRVparams.windowlength; % segment length 300s (5 min)

% try
    
   % Preprocess the PPG signal and extract IBI
   if ~exist([SavePath filesep 'FiltSigs' filesep SigName '_filt.mat'],'file')
       ppg_filt = get_filtered_PPG(PPG,HRVparams,1,SigName);       
   else
       load([SavePath filesep 'FiltSigs' filesep SigName '_filt.mat'],'ppg_filt');
   end
    
    % check if annotations are already presents
    if exist([SavePath filesep 'Annotation' filesep SigName '.sqippg'],'file')
        [onsets,~,sqinum] = read_ann([SavePath filesep 'Annotation' filesep SigName],'sqippg');
%%%%%        sqinum = sqinum(2:end); % convert sqi to range 0-100
    else
        % PPG Detection - qppg
        onsets = qppg_fast(ppg_filt,HRVparams.Fs);
        % PPG SQI 
        [sqinum, ~, ppgsqi] = calculate_ppgsqi(onsets,ppg_filt,HRVparams.Fs);
        % Write PPG  annotations
        AnnotationFolder = [SavePath filesep 'Annotation'];
        if ~exist(AnnotationFolder,'dir'); mkdir(AnnotationFolder); end
        write_ann([AnnotationFolder filesep SigName],HRVparams,'sqippg',onsets(1:length(sqinum)),char(ppgsqi),sqinum);
    end  
    
    RIAV = Extract_Respiratory_Components(ppg_filt,onsets,HRVparams.Fs); % respiratory induced  amplitude  variation
%%%%% adding sqi_resamp processing
    [NN_resamp, RIAV_resamp, sqi_resamp] = clean_resample_timeseries(onsets,RIAV,sqinum,HRVparams.Fs,Fs_resamp); 
    
    % segmnets of 5-min-length, sliding forwrard every 30s
    N_win = floor(length(RIAV_resamp)/Fs_resamp-win_len)/win_inc+1;
    
    
    for idx = 1:N_win
        % image data number
        
        try
            rr_in_win = NN_resamp((idx-1)*win_inc*Fs_resamp+1:min(((idx-1)*win_inc+win_len)*Fs_resamp,length(NN_resamp)));
            RIAV_in_win = RIAV_resamp((idx-1)*win_inc*Fs_resamp+1:min(((idx-1)*win_inc+win_len)*Fs_resamp,length(RIAV_resamp)));
            
            rr = rr_in_win;
            rrs = (1/Fs_resamp):(1/Fs_resamp):length(rr_in_win)/Fs_resamp;           

            %%% HRV parameters from NN intervals
            SDNN=std(rr);  
            
            % Using HRV_toolbox functions
            [AC,DC] = prsa(rr, rrs, HRVparams, [], 1);
            % SampEn
            sampEn_10 = EvalEntropyMetrics(rr, cumsum(rr),2,0.10, HRVparams, 1, []);
            sampEn_15 = EvalEntropyMetrics(rr, cumsum(rr),2,0.15, HRVparams, 1, []);
            sampEn_20 = EvalEntropyMetrics(rr, cumsum(rr),2,0.20, HRVparams, 1, []);
            % LF/HF
            [PSD,F] = CalcLomb(rrs',detrend(rr'), 1/NFFT:1/NFFT:(Fs_resamp/2) ,[],0);
            LF_HF_ratio=sum(PSD(F>0.04&F<=0.15))/sum(PSD(F>0.15&F<=0.4));

                 
            x = detrend(rr_in_win);
            y = detrend(RIAV_in_win);
            Pxy_amp = abs(csd_amp(x,y,2^10,Fs_resamp,hanning(512),256));
            Pxy_phase = abs(csd_phase(x,y,2^10,Fs_resamp,hanning(512),256));
            
            % get max amplitude and normalize
            max_amp = max(Pxy_amp);
            max_phase = max(Pxy_phase);
            z1 = (Pxy_amp./max_amp).^2;
            z2 = (Pxy_phase./max_phase).^2;
            Cxy = z1.*z2;
            
            % call crc % from Robert Joseph Thomas, An Electrocardiogram-Based
            % Technique to Assess Cardiopulmonary Coupling During Sleep
            
            cross_sp=cpsd(x,y,hanning(512),[],1024);
            co=coherence(x(1:1024),y(1:1024),hanning(512),1024);
            crc = abs(cross_sp).^2 .* co';
            crc=crc./max(crc);
            

            [low_peaks1, ~]=findpeaks(Cxy(low_band),'SORTSTR','descend');
            [high_peaks1, ~]=findpeaks(Cxy(high_band),'SORTSTR','descend');
            if ~isempty(low_peaks1) && ~isempty(low_peaks1)
                Lo_Hi_ratio1=sum(low_peaks1(1:min(2,length(low_peaks1))))/sum(high_peaks1(1:min(2,length(low_peaks1))));
            else
                if isempty(low_peaks1)
                    Lo_Hi_ratio1=0;
                end
                if isempty(high_peaks1)
                    Lo_Hi_ratio1=NaN;
                end
            end
            [low_peaks2, ~] = findpeaks(crc(low_band), 'SORTSTR' ,'descend');
            [high_peaks2, ~] = findpeaks(crc(high_band), 'SORTSTR' ,'descend');
            if ~isempty(low_peaks2) && ~isempty(low_peaks2)
                Lo_Hi_ratio2=sum(low_peaks2(1:min(2,length(low_peaks2))))/sum(high_peaks2(1:min(2,length(low_peaks2))));
            else
                if isempty(low_peaks2)
                    Lo_Hi_ratio2=0;
                end
                if isempty(high_peaks2)
                    Lo_Hi_ratio2=NaN;
                end
            end
            
            % area ratio
            sum_very_low=sum(crc(very_low_band));
            sum_low=sum(crc(low_band));
            sum_high=sum(crc(high_band));
            area_ratio_v_l=sum_very_low/sum_low;
            area_ratio_v_h=sum_very_low/sum_high;
            area_ratio_l_h=sum_low/sum_high;
            
            % set each image by a 512 window and slide forward every 10 seconds
            moving_window=512;
            sliding=10*Fs_resamp;
            for k=1:18
                x=rr_in_win((k-1)*sliding+1:(k-1)*sliding+moving_window);
                y=RIAV_in_win((k-1)*sliding+1:(k-1)*sliding+moving_window);
                x=detrend(x);
                y=detrend(y);
                
                % call crc % from Robert Joseph Thomas, An Electrocardiogram-Based
                % Technique to Assess Cardiopulmonary Coupling During Sleep
                
                cross_sp=cpsd(x,y,hanning(256),[],NFFT);
                co=coherence(x,y,hanning(256),NFFT);
                crc = abs(cross_sp).^2 .* co';
                crc = crc./max(crc);
                imagedata = crc(1:250);
                % resample crc data to 50 points
                tmpImg = [imagedata(1:25); resample(imagedata(26:end),25,225)];
                imdb.images.data(:,k,1,idx)= tmpImg;
            end
            
            % set class of image
            % SQI_i = sqinum((idx-1)*win_inc+1:min(((idx-1)*win_inc+win_len),length(sqinum)));
            SQI_i = nanmean(sqi_resamp((idx-1)*win_inc*Fs_resamp+1:min(((idx-1)*win_inc+win_len)*Fs_resamp,length(sqi_resamp))));
            
            imdb.meta.sqi(:,idx) = SQI_i; 
            imdb.meta.features(idx,1) = AC;
            imdb.meta.features(idx,2) = DC;
            imdb.meta.features(idx,3) = sampEn_10;
            imdb.meta.features(idx,4) = sampEn_15;
            imdb.meta.features(idx,5) = sampEn_20;
            imdb.meta.features(idx,6) = NaN;
            imdb.meta.features(idx,7) = SDNN;
            imdb.meta.features(idx,8) = LF_HF_ratio;
            imdb.meta.features(idx,9) = Lo_Hi_ratio1;
            imdb.meta.features(idx,10) = Lo_Hi_ratio2;
            imdb.meta.features(idx,11) = area_ratio_v_l;
            imdb.meta.features(idx,12) = area_ratio_v_h;
            imdb.meta.features(idx,13) = area_ratio_l_h;
            imdb.meta.rr{idx} = rr;
        catch
            
            error=1
            imdb.meta.sqi(:,idx) = NaN; 
            imdb.meta.features(idx,1) = NaN;
            imdb.meta.features(idx,2) = NaN;
            imdb.meta.features(idx,3) = NaN;
            imdb.meta.features(idx,4) = NaN;
            imdb.meta.features(idx,5) = NaN;
            imdb.meta.features(idx,6) = NaN;
            imdb.meta.features(idx,7) = NaN;
            imdb.meta.features(idx,8) = NaN;
            imdb.meta.features(idx,9) = NaN;
            imdb.meta.features(idx,10) = NaN;
            imdb.meta.features(idx,11) = NaN;
            imdb.meta.features(idx,12) = NaN;
            imdb.meta.features(idx,13) = NaN;
%             imagedata = zeros(1,250);
            imdb.meta.rr{idx} = zeros(1,1200)*NaN;%[imagedata(1:25); resample(imagedata(26:end),25,225)];
            continue;
        end
    end
        
% catch
% end
imdb.images.data=single(imdb.images.data);

end


