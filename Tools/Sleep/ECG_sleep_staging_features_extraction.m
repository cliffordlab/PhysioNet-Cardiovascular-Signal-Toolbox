function features = ECG_sleep_staging_features_extraction(ECGdata,fs)

% try
    NFFT=2000;
    imdb_n=0;

    ecg_data=ECGdata;
    Fs=fs;
    
    ecg_data(find(isnan(ecg_data)))=nanmean(ecg_data);
    
    clear ECGdata
    
    if Fs~=125
        ecg_data=resample(ecg_data,125,Fs);
        Fs=125;
    end
    
    % filter comparison
    ecg_filter1=ecg_LP_filter_125_22(ecg_data);
    ecg_filter1=ecg_HP_filter_125_012(ecg_filter1);
    
    clear ecg_data
    % set to positive peaks
    pos_p=[];
    neg_p=[];
    seg_p=1;
    for j=1:Fs*5:length(ecg_filter1) % find peaks every 5 sec
        pos_p(seg_p)=max(ecg_filter1(j:min(length(ecg_filter1),j+Fs*5-1)));
        neg_p(seg_p)=min(ecg_filter1(j:min(length(ecg_filter1),j+Fs*5-1)));
        seg_p=seg_p+1;
    end
    if abs(median(pos_p))<abs(median(neg_p))
        ecg_filter1=-ecg_filter1;
    end
    
    % if the amplitude of QRS is too small, enlarge it
    QRS_amp=median(pos_p)-median(neg_p);
    if QRS_amp<1.0
        ecg_filter1=ecg_filter1.*(1.0/QRS_amp);
    end
    
    % QRS detection
    qrs_pos1 = run_qrsdet_by_seg3(ecg_filter1,Fs,15,0.6,'MECG');
    
    for (j=1:length(qrs_pos1)) % for 125 Hz, find +/- 7 points
        [pk, ploc]=findpeaks(ecg_filter1(max(1,qrs_pos1(j)-7+1):min(length(ecg_filter1),qrs_pos1(j)+7)));
        if length(pk)>0
            [maxpk,maxloc]=max(pk);
            qrs_pos1(j)=max(1,qrs_pos1(j)-7+1)+ploc(maxloc)-1;
        end
    end
    
    % QRS detection by wqrs
    qrs_pos2 = wqrsm_fast(ecg_filter1,Fs);
    
    % bsqi
    bsqi=zeros(1,floor(length(ecg_filter1)/Fs));
    ann1=qrs_pos1/Fs;
    ann2=qrs_pos2/Fs;
    endtime=max([ann1(end) ann2(end)]);
    for j=0:round(endtime)-10
        [sqi]=run_sqi_n(ann1,ann2,0.1,j,10,Fs);
        if j==0
            bsqi(1:j+5)=round(sqi*100);
        else
            bsqi(j+5)=round(sqi*100);
        end
    end
    bsqi(j+5:j+10)=round(sqi*100);
    bsqi(find(isnan(bsqi)))=0;
    SQI=bsqi;
    
    qrs_pos=qrs_pos1;
    
    % get the R-R interval and QRS amplitude for each beat
    RR1=diff(qrs_pos1)./Fs;
    QRS_A1=ecg_filter1(qrs_pos1(2:end));
    %     y_edr1=edr(0,ecg_filter,qrs_pos/Fs,Fs);
    %     QRS_A=y_edr(length(y_edr)/2+2:end);
    qrs_pos1=qrs_pos1(2:end);
    
    %     % reject ectopic beats, was changed to new algorithm below
    %     RR_e1 = ectopic_rejection(RR1);
    %     valid_RR1=find(RR_e1~=-1);
    %     qrs_pos1=qrs_pos1(valid_RR1);
    %     QRS_A1=QRS_A1(valid_RR1);
    %     RR_e1=RR_e1(valid_RR1);
    
    clear ecg_filter1
    
    RR_e0=RR1;    % no ectopic rejection
    QRS_A0=QRS_A1;
    
    % reject ectopic beats
    % remove invalid R-R intervals if RR<0.3 or RR>2.0 & 20% outside of 41 points moving mean exclude
    RR_ori = RR1;
    RR_mean_41=meanfilt1(RR_ori,41);
    RR=[]; % RR interval in seconds
    RR_time=[]; % time of RR interval in seconds
    v_RR=ones(1,length(RR_ori))*-1;
    n_RR=0;
    for j=1:length(RR_ori)
        if abs(RR_ori(j)-RR_mean_41(j)) < 0.2*RR_mean_41(j) & ...
                RR_ori(j)>= 0.3 & RR_ori(j) <= 2.0
            n_RR=n_RR+1;
            RR(n_RR)=RR_ori(j);
            RR_time(n_RR)=qrs_pos1(j)/Fs;
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
            QRS(n_QRS)=QRS_A1(k);
            QRS_time(n_QRS)=qrs_pos1(k)/Fs;
            v_QRS(k)=1;
        end
    end
    valid_QRS=find(v_QRS==1);
    valid_RRQRS=intersect(valid_RR,valid_QRS);
    valid_RR=valid_RRQRS;
    
    qrs_pos1=qrs_pos1(valid_RR);
    QRS_A1=QRS_A1(valid_RR);
    RR_e1=RR1(valid_RR);
    
    % resample to 2Hz ??? 20Hz
    qrs_pos1=qrs_pos1./Fs;
    t2=round(qrs_pos1(end));
    Fs_resamp=4;
    t1=(1/Fs_resamp):(1/Fs_resamp):t2; % 4Hz resample
    
    QRS_A_resamp=interp1(qrs_pos1,QRS_A1,t1);%,'spline');
    RR_e_resamp=interp1(qrs_pos1,RR_e1,t1);%,'spline');
    QRS_A_resamp(find(isnan(QRS_A_resamp)))=nanmean(QRS_A_resamp);
    RR_e_resamp(find(isnan(RR_e_resamp)))=nanmean(RR_e_resamp);
    
    SQI_resamp=interp1(1:length(SQI),SQI,1:1/Fs_resamp:length(SQI));
    SQI_resamp=[SQI(1) SQI(1) SQI(1) SQI_resamp];
      
%     % stage_n: 1:wake, 2:rem, 3:s1, 4:s2, 5:s3
%     load([annfolder filesep filename '-ann.mat']);
%     % stages previous type: 0 - wake, 1~4 - sleep 1~4, 5 - REM
%     stage=stages;
%     stage=stage+1;
%     % change to stage_n
%     stage(find(stage==2))=3;
%     stage(find(stage==5))=4;
%     stage(find(stage==6))=2;
%     % after: 1:wake, 2:rem, 3:nrem light, 4:nrem deep
%     
%     stage_n=stage;
%     stage_time=[0:length(stage_n)-1]*30;
%     stage_type=[4 3 2 1];
    
    ann_epoch_len=30; % annotation epoch is 30s
    seg_length=300; % segment length 300s (5 min)
    epochs_per_seg=seg_length/ann_epoch_len;
    slide_window=30;
    
%     clear stage
    HRVparams = InitializeHRVparams([]);
    % Using HRV_toolbox functions
    % PRSA
    HRVparams.prsa.thresh_per=5;
    HRVparams.prsa.win_length=2;
    HRVparams.plot_results=0;
    HRVparams.prsa.scale=2;
    HRVparams.prsa.plot_results=0;
    HRVparams.windowlength=300;
    
    HRVparams.sqi.LowQualityThreshold=0.5;
    HRVparams.RejectionThreshold=0.5;
    HRVparams.af.on=0;

    HRVparams.output.separate = 0;
    HRVparams.MSE.on = 0;
    HRVparams.DFA.on = 0;
    HRVparams.HRT.on =0;
    HRVparams.output.format='mat';
    
    % HRV parameters
    
    % duplicate the first 30s-segment and last 30-s segment 5 times (150s) for padding
    QRS_A_resamp=[QRS_A_resamp(1:120) QRS_A_resamp(1:120)  QRS_A_resamp(1:120) QRS_A_resamp(1:120) QRS_A_resamp(1:120) QRS_A_resamp QRS_A_resamp(end-120+1:end) QRS_A_resamp(end-120+1:end) QRS_A_resamp(end-120+1:end)  QRS_A_resamp(end-120+1:end) QRS_A_resamp(end-120+1:end)];
    RR_e_resamp=[RR_e_resamp(1:120) RR_e_resamp(1:120)  RR_e_resamp(1:120) RR_e_resamp(1:120) RR_e_resamp(1:120) RR_e_resamp RR_e_resamp(end-120+1:end) RR_e_resamp(end-120+1:end) RR_e_resamp(end-120+1:end)  RR_e_resamp(end-120+1:end) RR_e_resamp(end-120+1:end)];
    SQI_resamp=[SQI_resamp(1:120) SQI_resamp(1:120)  SQI_resamp(1:120) SQI_resamp(1:120) SQI_resamp(1:120) SQI_resamp SQI_resamp(end-120+1:end) SQI_resamp(end-120+1:end) SQI_resamp(end-120+1:end)  SQI_resamp(end-120+1:end) SQI_resamp(end-120+1:end)];
    t2=t2+10*30;
    t1=(1/Fs_resamp):(1/Fs_resamp):t2; % 4Hz resample
  
    
    seg_select=6;
    stage_n=zeros(1,floor(length(QRS_A_resamp)/ann_epoch_len/Fs_resamp));
    
    rr=RR_e_resamp((seg_select-1)*ann_epoch_len*Fs_resamp-135*Fs_resamp+1:end);
    t_rr = t1((seg_select-1)*ann_epoch_len*Fs_resamp-135*Fs_resamp+1:end);
    sqi_rr=SQI_resamp((seg_select-1)*ann_epoch_len*Fs_resamp-135*Fs_resamp+1:end);
    
    sqi_len=min([length(rr) length(t_rr) length(sqi_rr)]);
    rr=rr(1:sqi_len);
    t_rr=t_rr(1:sqi_len);
    sqi_rr=sqi_rr(1:sqi_len);
    
    sqi_rr=[t_rr;sqi_rr];
    
    if size(sqi_rr,1)<size(sqi_rr,2)
        sqi_rr=sqi_rr';
    end
    [HRVout] = Main_HRV_Analysis(rr,t_rr,'RRIntervals',HRVparams,' ',[],sqi_rr);

    
    for seg_select=6:length(stage_n)-5
        %         if ~isempty(intersect(stage_type,stage_n(seg_select))) % valid sleep stage annotation
        % image data number
                       
        imdb_n=imdb_n+1;

        try
            % select 5-min epochs by centered the 30s of annotation
                seg=(seg_select-1)*ann_epoch_len*Fs_resamp-135*Fs_resamp+1:(seg_select)*ann_epoch_len*Fs_resamp+135*Fs_resamp;
                
                x_seg=RR_e_resamp(seg);
                y_seg=QRS_A_resamp(seg);
                
                rr=x_seg;
                rrs=(1/Fs_resamp):(1/Fs_resamp):length(x_seg)/Fs_resamp;
                % SDNN
                SDNN=std(rr);
                
                % Using HRV_toolbox functions
                % PRSA
                %                 HRVparams.prsa.thresh_per=5;
                %                 HRVparams.prsa.win_length=2;
                %                 HRVparams.plot_results=0;
                %                 HRVparams.prsa.scale=2;
                %                 HRVparams.prsa.plot_results=0;
                %                 HRVparams.windowlength=300;
                %                 HRVparams.sqi.LowQualityThreshold=0.8;
                %                 HRVparams.RejectionThreshold=0.2;
                [AC,DC, prsa_ac, prsa_dc] = prsa(rr, rrs, HRVparams, [], 1);
                % MSE
                y10 = ComputeMultiscaleEntropy(rr',2,0.10,1);
                y15 = ComputeMultiscaleEntropy(rr',2,0.15,1);
                y20 = ComputeMultiscaleEntropy(rr',2,0.20,1);
                %             % DFA
                %             DFA_alpha=dfaScalingExponent(rr);
                % LF/HF
                [PSD,F] = CalcLomb(rrs',detrend(rr'),[1/NFFT:1/NFFT:(Fs_resamp/2)],[],0);
                LF_HF_ratio=sum(PSD(find(F>0.04&F<=0.15)))/sum(PSD(find(F>0.15&F<=0.4)));
                
                x=detrend(x_seg);
                y=detrend(y_seg);
                % call cpc % from Laurence Zapanta, Heart Rate Variability in Mice with Coronary Heart Disease
                % [CRP,freq] = cohere_csd_lz(RR_e_2Hz,QRS_A_2Hz,2^9,Fs,hanning(400),200)
                [Pxy_amp, Fxy] = csd_amp(x,y,2^10,Fs_resamp,hanning(512),256);
                [Pxy_phase, Fxy] = csd_phase(x,y,2^10,Fs_resamp,hanning(512),256);
                
                Pxy_phase = abs(Pxy_phase);
                Pxy_amp = abs(Pxy_amp);
                
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
                
                % 512 points for 2 Hz (Fs_resamp/2)
                very_low_band=1:3;  % <0.01
                low_band=4:25;      % 0.01-0.1
                high_band=26:103;   % 0.1-0.4
                [low_peaks1 low_locs1]=findpeaks(Cxy(low_band),'SORTSTR','descend');
                [high_peaks1 high_locs1]=findpeaks(Cxy(high_band),'SORTSTR','descend');
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
                [low_peaks2 low_locs2]=findpeaks(crc(low_band),'SORTSTR','descend');
                [high_peaks2 high_locs2]=findpeaks(crc(high_band),'SORTSTR','descend');
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
                
                imagedata=crc(1:128); % 0-0.5Hz
                
%                 % HRV parameters
%                 t_rr = rrs;
%                 sqi_rr=SQI_resamp(seg);
% 
%                 sqi_rr=[t_rr;sqi_rr];
%                 
%                 if size(sqi_rr,1)<size(sqi_rr,2)
%                     sqi_rr=sqi_rr';
%                 end
%                 [HRVout] = Main_HRV_Analysis(rr(1:end-1),t_rr(1:end-1),'RRIntervals',HRVparams,'',[],sqi_rr(1:end-1,:));
                
                 
                % resample crc data to 50 points
                imdb.images.data(:,imdb_n)=[imagedata(1:25); resample(imagedata(26:end),25,103)];
                
                imdb.images.labels(imdb_n)=stage_n(seg_select);
                imdb.meta.sqi(:,imdb_n)=SQI((seg_select-1)*ann_epoch_len-135+1:(seg_select)*ann_epoch_len+135);
                
                imdb.meta.features(imdb_n,1)=AC;
                imdb.meta.features(imdb_n,2)=DC;
                imdb.meta.features(imdb_n,3)=y10;
                imdb.meta.features(imdb_n,4)=y15;
                imdb.meta.features(imdb_n,5)=y20;
                imdb.meta.features(imdb_n,6)=SDNN;
                imdb.meta.features(imdb_n,7)=LF_HF_ratio;
                imdb.meta.features(imdb_n,8)=Lo_Hi_ratio1;
                imdb.meta.features(imdb_n,9)=Lo_Hi_ratio2;
                imdb.meta.features(imdb_n,10)=area_ratio_v_l;
                imdb.meta.features(imdb_n,11)=area_ratio_v_h;
                imdb.meta.features(imdb_n,12)=area_ratio_l_h;
                imdb.meta.rr{imdb_n}=rr;
                imdb.meta.rrs{imdb_n}=rrs;
%                 imdb.meta.HRV_features(:,imdb_n)=HRVout;
        catch
                imdb.images.data(:,imdb_n)=zeros(50,1).*NaN;
                
                imdb.images.labels(imdb_n)=stage_n(seg_select);
                imdb.meta.sqi(:,imdb_n)=zeros(300,1).*NaN;
                
                imdb.meta.features(imdb_n,1)=NaN;
                imdb.meta.features(imdb_n,2)=NaN;
                imdb.meta.features(imdb_n,3)=NaN;
                imdb.meta.features(imdb_n,4)=NaN;
                imdb.meta.features(imdb_n,5)=NaN;
                imdb.meta.features(imdb_n,6)=NaN;
                imdb.meta.features(imdb_n,7)=NaN;
                imdb.meta.features(imdb_n,8)=NaN;
                imdb.meta.features(imdb_n,9)=NaN;
                imdb.meta.features(imdb_n,10)=NaN;
                imdb.meta.features(imdb_n,11)=NaN;
                imdb.meta.features(imdb_n,12)=NaN;
                imdb.meta.rr{imdb_n}=NaN;
                imdb.meta.rrs{imdb_n}=NaN;
%                 imdb.meta.HRV_features(:,imdb_n)=zeros(29,1).*NaN;
            
        end

    end
    imdb_n
    while size(HRVout,1)<imdb_n
        HRVout=[HRVout;zeros(1,29).*NaN];
    end
    imdb.meta.HRV_features=HRVout(1:imdb_n,:);
    
    nfeatures_features=12;
    nfeatures_HRV_used=[3:5 7:12 14 16:21 23:29];
    
    for l=7:nfeatures_features
        for j=1:size(imdb.meta.features,1)
            imdb.meta.features(j,l)=log(imdb.meta.features(j,l));
        end
    end
    imdb.meta.features(find(isinf(imdb.meta.features)))=NaN;
    
    sqi=mean(imdb.meta.sqi,1);
        
    features=[imdb.images.data;imdb.meta.features';sqi;imdb.meta.HRV_features(:,nfeatures_HRV_used)'];
    
% catch
%     display(['Feature extraction Error!!!']);
%     features=[];
% end
end
