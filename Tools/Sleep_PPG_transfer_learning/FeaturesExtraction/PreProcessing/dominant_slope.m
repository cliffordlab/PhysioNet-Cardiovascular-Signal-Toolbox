function domi_factor = dominant_slope(data,fs)
% decide whch side has the dominant slope (up or down)
%
% return domi_factor = 1, up slope is dominant
%                     -1, down slope is dominant
% so the data may be multiplied by domi_factor

domi_factor=1;
up_slope=0;
down_slope=0;
% split data to every 2s and find the median amplitude of each segment
d=reshape(data(1:floor(length(data)/(fs*2))*fs*2),fs*2,[]);
amplitude=median(max(d)-min(d));
% use amplitude/3 as MinPeakProminence to find peaks
[pks locs]=findpeaks(data,'MinPeakProminence',amplitude/3);
% calculate x-coordinate of each slope
for i=1:length(pks)-1
    [mn mni]=min(data(locs(i):locs(i+1)));
    l=length(data(locs(i):locs(i+1)));
    up_slope(i)=l-mni;
    down_slope(i)=mni-1;
end
% sum(up_slope)
% sum(down_slope)
% median(up_slope)
% median(down_slope)
if median(up_slope)>median(down_slope)
    domi_factor=-1;
end
if median(up_slope)==median(down_slope)
    if sum(up_slope)>sum(down_slope)
        domi_factor=-1;
    end
end
end
    
    