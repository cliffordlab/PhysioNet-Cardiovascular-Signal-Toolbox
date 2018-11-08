% readdata.m read PPG data and annotation from database record from the 
% begintime to the stoptime 

function [d ann a]=readdata(datanumber,annname,PPG_lead,begintime,stoptime)
if nargin<5
    stoptime = [];
end
if nargin<4
    begintime = ['00:00:00'];
end
if nargin<3
    fprintf('Error: must provide datanumber, annotation name and PPG_lead');
    return;
end

if length(stoptime)==0
    r=rdsamp(datanumber,'begin',begintime,'sigs',PPG_lead);
else
    r=rdsamp(datanumber,'begin',begintime,'stop',stoptime,'sigs',PPG_lead);
end
    
if length(r)>0
    d=r(:,2);
    for i=1:length(d)
        if d(i)<-32700
            if (i)==1
                d(i)=0;
            else
                d(i)=d(i-1);
            end
        end
    end
else
    d=[];
end
if length(stoptime)==0
    a=rdann(datanumber,annname,'start',begintime);
else
    a=rdann(datanumber,annname,'start',begintime,'stop',stoptime);
end
if length(a)<1
    ann=[];
end
for i=1:length(a)
    ann(i)=a(i).sampleNumber;
end
