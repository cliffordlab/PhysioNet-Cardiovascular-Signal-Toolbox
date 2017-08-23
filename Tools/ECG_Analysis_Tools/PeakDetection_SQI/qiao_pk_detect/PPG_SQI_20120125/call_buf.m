function call_buf()
% Make sure you can access the MIMIC database when run the example
% It will read the 039/039 waveform data and 039.wpleth or 039.ple annotation
% and create the PPGSQI annotation named 039.wplesqi
% Make sure you have 'The WFDB Toolbox for Matlab' installed
% see: http://www.physionet.org/physiotools/matlab/wfdb-swig-matlab/

filename='039/039';
annname='wpleth'; % annname='ple';
PPG_lead=2; % lead number of Pleth in 039
b='0:0:0'; % begin time to analyse
s='0:10:0'; % end time to analyse
samp_freq=125;

% Read data from WFDB data file
d=rdsamp(filename,'begin',b,'stop',s,'sigs',PPG_lead);
d=d(:,2);

% PPG beat detection
[beat] = wabp_pleth_new(d);
        
if length(d)<1 || length(beat)<1
     fprintf('Error: read PPG data or PPG beat detection error\n');   
     return;
end

% analysis PPG_SQI          
windowlen=30*samp_freq; % 30s window             
template=[];           
beat_i=1;
            
for j=1:ceil(length(d)/windowlen)                
    databegin=(j-1)*windowlen+1;                
    dataend=min(length(d),j*windowlen);                
    annf=find(beat<=dataend);                
    if length(annf)<1                    
        continue;                
    end    
    if length(annf)<length(beat)                    
        annf(length(annf)+1)=annf(length(annf))+1;                
    end    
    annf=find(beat(annf)>=databegin);                
    if length(annf)<1                    
        continue;                
    end    
    anntime=beat(annf);                
    wave=d(databegin:min(length(d),max(dataend,anntime(length(anntime))+3*samp_freq))); % prolong the data window an extra 3s                
    % set the anntime to be the 0 offset of the selected data
    anntime=beat(annf)-databegin+1;         
    % PPG SQI analysis
    [annot sqimatrix template valid] = PPG_SQI_buf(wave,anntime,template,30*samp_freq);                
    for i=1:length(annot)                    
        annot_all(beat_i)=annot(i);                    
        sqimatrix_all(beat_i,:)=sqimatrix(i,:);                    
        beat_i=beat_i+1;                
    end
end

analysisperiod_s=floor(length(d)/samp_freq); % data length (s)
moving_window=10; % moving window (s) for SQI analysis

% transfer PPG SQI marker to 1/0
qsqi=[];
if length(annot_all)>0
        for k=1:length(annot_all)
            if annot_all(k)=='E' || annot_all(k)=='A'
                qsqi(k)=1;
            else
                qsqi(k)=0;
            end
        end
end

% Second-by-second PPG SQI
for k=1:analysisperiod_s-moving_window+1
    if k==1
        sbs_qsqi(1:moving_window)=find_mean_sqi(beat,qsqi,0,moving_window);
        l=moving_window+1;
    else
        sbs_qsqi(l)=find_mean_sqi(beat,qsqi,k-1,k+moving_window-1);
        l=l+1;
    end
end

figure
hold on
plot(d);
plot(beat(1:length(qsqi)),qsqi*1000,'g');
plot(beat,0,'rx')
legend('PPG wave','PPG SQI','Beat detection')

                   
function sqi_out = find_mean_sqi(ann, sqi, begin_s, end_s, samp_freq)
% calculate the mean sqi within begin_s to end_s seconds (jsqi,qsqi)
%
% input:
%   ann       : annotation sampleNumber offset of each beat
%   sqi       : jsqi/qsqi buffur
%   begin_s   : begin of calculate window of the annotation (s)
%   end_s     : end of calculate window of the annotation (s)
%   samp_freq : sampling frequency, default : 125 Hz
%
% output:
%   sqi_out   : mean sqi value within the window

sqi_out=0;

if (nargin < 5)
    samp_freq = 125;
end

begin_samp = begin_s*samp_freq;
end_samp = end_s*samp_freq;

first_ann = find(ann>=begin_samp,1,'first');
last_ann =  find(ann<=end_samp,1,'last');
if last_ann - first_ann > 0
    sqi_out = mean(sqi(first_ann:min(length(sqi),last_ann)));
end

