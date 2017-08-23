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


        
	[d ann a]=readdata(filename,annname,PPG_lead,b,s);
        if length(d)<1 || length(ann)<1
            fprintf('Error: read PPG data or wpleth annotation error\n');
        end
            % analysis PPG_SQI
            windowlen=30*samp_freq; % 30s window 
            template=[];
%             ann=ann-(alarmInd(i)-(60+13)*samp_freq);
            for j=1:ceil(length(d)/windowlen)
                databegin=(j-1)*windowlen+1;
                dataend=min(length(d),j*windowlen);
                annf=find(ann<=dataend);
                if length(annf)<1
                    continue;
                end
                if length(annf)<length(ann)
                    annf(length(annf)+1)=annf(length(annf))+1;
                end
                annf=find(ann(annf)>=databegin);
                if length(annf)<1
                    continue;
                end
                anntime=ann(annf);
                annot=a(annf);
                wave=d(databegin:min(length(d),max(dataend,anntime(length(anntime))+3*samp_freq))); % prolong the data window an extra 3s
                anntime=ann(annf)-databegin+1;
                [annot template valid] = PPG_SQI(wave,anntime,annot,template,30*samp_freq);
                a(annf)=annot;
            end
            

            wrann(a,filename,'wplesqi');

