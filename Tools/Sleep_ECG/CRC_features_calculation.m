function imdb = CRC_features_calculation(ECG,OP)

imdb=[];
imdb.images.data=[];
imdb.mata.feature=[];
try
    NFFT=2000;
    imdb_n=0;
    Fs=OP.fs;
    
    ECG(find(isnan(ECG)))=nanmean(ECG);
    
    % preprecessing
    ecg_filter1=OP.lp_filter(ECG);
    ecg_filter1=OP.hp_filter(ecg_filter1);

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
    
    % QRS detection
	qrs_pos1 = run_qrsdet_by_seg3(ecg_filter1,Fs,15,0.6,'MECG');
    
    % search span for max of QRS peak between qrs_pos1 +/- qrs_span
    qrs_span=round(Fs*0.05);
    
    for (j=1:length(qrs_pos1))
        [pk, ploc]=findpeaks(ecg_filter1(max(1,qrs_pos1(j)-qrs_span+1):min(length(ecg_filter1),qrs_pos1(j)+qrs_span)));
        if length(pk)>0
            [maxpk,maxloc]=max(pk);
            qrs_pos1(j)=max(1,qrs_pos1(j)-qrs_span+1)+ploc(maxloc)-1;
        end
    end
    
    % QRS detection by wqrs
    qrs_pos2 = wqrsm(ecg_filter1,Fs);
    
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
    qrs_pos1=qrs_pos1(2:end);
      
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
    
    % resample to 4Hz
    qrs_pos1=qrs_pos1./Fs;
    t2=round(qrs_pos1(end));
    Fs_resamp=4;
    t1=(1/Fs_resamp):(1/Fs_resamp):t2; % 4Hz resample
      
    % linear interpolation
    QRS_A_resamp=interp1(qrs_pos1,QRS_A1,t1);%,'spline');
    RR_e_resamp=interp1(qrs_pos1,RR_e1,t1);%,'spline');
    QRS_A_resamp(find(isnan(QRS_A_resamp)))=nanmean(QRS_A_resamp);
    RR_e_resamp(find(isnan(RR_e_resamp)))=nanmean(RR_e_resamp);


    sliding_win_len=30; % sliding window length is 30s
    seg_length=300; % segment length 300s (5 min)
    
    % segmnets of 5-min-length, sliding forwrard every 30s
    n_seg=floor(length(QRS_A_resamp)/Fs_resamp-seg_length)/sliding_win_len+1;
    
    for j=1:n_seg
        % image data number
        try
            imdb_n=imdb_n+1;
            
            x_seg=RR_e_resamp((j-1)*sliding_win_len*Fs_resamp+1:min(((j-1)*sliding_win_len+seg_length)*Fs_resamp,length(RR_e_resamp)));
            y_seg=QRS_A_resamp((j-1)*sliding_win_len*Fs_resamp+1:min(((j-1)*sliding_win_len+seg_length)*Fs_resamp,length(RR_e_resamp)));
            
            rr=x_seg;
            rrs=(1/Fs_resamp):(1/Fs_resamp):length(x_seg)/Fs_resamp;
%            HR_seg=60./x_seg;

            % SDNN
            SDNN=std(rr);
            
%             % PRSA
%             [acm0 dcm0 AC DC] = prsa_gari(rr,rrs,[],2,5,0);
%             % MSE
%             y10=msentropy(rr',[],[],[],[],[],[],[],1,0.1,0.1);
%             y15=msentropy(rr',[],[],[],[],[],[],[],1,0.15,0.15);
%             y20=msentropy(rr',[],[],[],[],[],[],[],1,0.2,0.2);
% %             % DFA
% %             [ln, lf]=dfa(rr');
% %             ft_=fittype('poly1');
% %             cf_=fit(ln,lf,ft_);
% %             DFA_alpha=cf_.p1;
%             % LF/HF
%             [Pxx,F]=lomb([rrs' rr']);
%             LF_HF_ratio=sum(Pxx(find(F>0.04&F<=0.15)))/sum(Pxx(find(F>0.15&F<=0.4)));
%
% Using HRV_toolbox functions
            % PRSA
            HRVparams.prsa.thresh_per=5;
            HRVparams.prsa.win_length=2;
            HRVparams.plot_results=0;
            HRVparams.prsa.scale=2;
            HRVparams.prsa.plot_results=0;
            HRVparams.windowlength=300;
            HRVparams.sqi.LowQualityThreshold=0.8;
            HRVparams.RejectionThreshold=0.2;
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
            
            % set each image by a 512 window and slide forward every 10 seconds
            moving_window=512;
            sliding=10*Fs_resamp;
            for k=1:18
                x=x_seg((k-1)*sliding+1:(k-1)*sliding+moving_window);
                y=y_seg((k-1)*sliding+1:(k-1)*sliding+moving_window);
                x=detrend(x);
                y=detrend(y);

%                 % call cpc % from Laurence Zapanta, Heart Rate Variability in Mice with Coronary Heart Disease
%                 % [CRP,freq] = cohere_csd_lz(RR_e_resamp,QRS_A_resamp,2^9,Fs,hanning(400),200)
%                 [Pxy_amp, Fxy] = csd_amp(x,y,NFFT,Fs_resamp,hanning(256),128);
%                 [Pxy_phase, Fxy] = csd_phase(x,y,NFFT,Fs_resamp,hanning(256),128);
%                 
%                 Pxy_phase = abs(Pxy_phase);
%                 Pxy_amp = abs(Pxy_amp);
%                 
%                 % get max amplitude and normalize
%                 max_amp = max(Pxy_amp);
%                 max_phase = max(Pxy_phase);
%                 z1 = (Pxy_amp./max_amp).^2;
%                 z2 = (Pxy_phase./max_phase).^2;
%                 Cxy = z1.*z2;
                
                % call crc % from Robert Joseph Thomas, An Electrocardiogram-Based
                % Technique to Assess Cardiopulmonary Coupling During Sleep
                
                cross_sp=cpsd(x,y,hanning(256),[],NFFT);
                co=coherence(x,y,hanning(256),NFFT);
                crc = abs(cross_sp).^2 .* co';
                crc=crc./max(crc);
                imagedata=crc(1:250);
                % resample crc data to 50 points
                imdb.images.data(:,k,1,imdb_n)=[imagedata(1:25); resample(imagedata(26:end),25,225)];
            end
            
            % set class of image
            SQI_i=SQI((j-1)*sliding_win_len+1:min(((j-1)*sliding_win_len+seg_length),length(SQI)));
            imdb.meta.sqi(:,imdb_n)=SQI_i; %((seg_begin(j)-1)*ann_epoch_len+1:(seg_begin(j)-1)*ann_epoch_len+seg_length);
            
            imdb.meta.features(imdb_n,1)=AC;
            imdb.meta.features(imdb_n,2)=DC;
            imdb.meta.features(imdb_n,3)=y10;
            imdb.meta.features(imdb_n,4)=y15;
            imdb.meta.features(imdb_n,5)=y20;
            imdb.meta.features(imdb_n,6)=NaN;
            imdb.meta.features(imdb_n,7)=SDNN;
            imdb.meta.features(imdb_n,8)=LF_HF_ratio;
            imdb.meta.features(imdb_n,9)=Lo_Hi_ratio1;
            imdb.meta.features(imdb_n,10)=Lo_Hi_ratio2;
            imdb.meta.features(imdb_n,11)=area_ratio_v_l;
            imdb.meta.features(imdb_n,12)=area_ratio_v_h;
            imdb.meta.features(imdb_n,13)=area_ratio_l_h;
            imdb.meta.rr{imdb_n}=rr;
        catch
            continue;
        end
    end
catch
end
imdb.images.data=single(imdb.images.data);

end


