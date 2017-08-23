%ectopic deletion
%delete next beat if the next beat is greater than the 
% previous beat by a threshold percentage of the previous beat

function [t,rr_out, ramp_out] = ectopic_rejection(t_hrv, rr, ramp, th)

len = length(rr);
i = 1;
j = 1;
t = t_hrv;
try
while i < len+1

    if rr(i+1) > (1+(j*th))*rr(i) | rr(i+1) < (1-(j*th))*rr(i)

        rr(i+1:i+2) = [];
        t(i+1:i+2) = [];
        ramp(i+1:i+2) = [];
        
        j = j+1;
        i = i-1;
    else
        j = 1;
    end
    i = i+1;
end
catch
end

ramp_out = ramp;
rr_out = rr;
n = len-length(rr_out);