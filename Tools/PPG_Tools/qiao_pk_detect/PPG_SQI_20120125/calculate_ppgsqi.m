function [ppgsqi_numeric sqi_mean_30 annot_all]=calculate_ppgsqi(PPGann,Waveform,fs)
        sqimatrix_all=zeros(length(PPGann),3);
        sqi_mean_30=zeros(1,ceil(length(Waveform)/fs/30));
        annot_all=[];
        for j=1:length(PPGann)
            annot_all(j)='Q';
        end
        
        % analysis PPG_SQI
        windowlen=30*fs; % 30s window
        template=[];
        beat_i=1;
        
        % loop every 30-sec
        for j=1:ceil(length(Waveform)/windowlen)
            try
                databegin=(j-1)*windowlen+1;
                dataend=min(length(Waveform),j*windowlen);
                annf=find(PPGann<=dataend);
                if length(annf)<=1
                    continue;
                end
                if length(annf)<length(PPGann)
                    annf(length(annf)+1)=annf(length(annf))+1;
                end
                annf=find(PPGann(annf)>=databegin);
                if length(annf)<=1
                    continue;
                end
                anntime=PPGann(annf);
                wave=Waveform(databegin:min(length(Waveform),max(dataend,anntime(length(anntime))+3*fs))); % prolong the data window an extra 3s
                % set the anntime to be the 0 offset of the selected data
                anntime=PPGann(annf)-databegin+1;
                % PPG SQI analysis
                [annot sqimatrix template valid] = PPG_SQI_buf(wave,anntime,template,30*fs,fs);
                for k=1:length(annot)
                    annot_all(annf(k))=annot{k};
                    sqimatrix_all(annf(k),:)=sqimatrix(k,1:3); % 1:4
                    beat_i=beat_i+1;
                end
                sqi_mean_30(j)=mean(mean(sqimatrix(:,1:3)'));
            catch
                fprintf('PPGsqi error at segment %d ',j);
            end
        end
        ppgsqi_numeric = round(mean(sqimatrix_all(:,1:3),2)');
end
